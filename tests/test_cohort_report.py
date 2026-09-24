#!/usr/bin/env python3
"""Checks for scripts/build_cohort_report.py, on synthetic rows only.

The report is the artefact most likely to be looked at and least likely to be tested,
so the checks here are the ones that would embarrass us: a page that needs the network
to render, a chart that overflows its own viewBox, a run on a dirty tree that does not
say so, and an empty or partial cohort that crashes instead of degrading.
"""

import csv
import importlib.util
import json
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_cohort_report", ROOT / "scripts" / "build_cohort_report.py"
)
CR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CR)

QC_SPEC = importlib.util.spec_from_file_location(
    "build_qc_summary", ROOT / "scripts" / "build_qc_summary.py"
)
QC = importlib.util.module_from_spec(QC_SPEC)
QC_SPEC.loader.exec_module(QC)


def sample_row(name="S1", **over):
    row = {c: "NA" for c in QC.COLUMNS}
    row.update(
        sample=name, reads_raw="1000000", reads_after_fastp="900000",
        pct_pass_fastp="90.00", duplication_rate="0.12",
        adapter_scan_applied="false",
        reads_after_hg38="200000", reads_after_t2t_phix="150000",
        reads_host_depleted="100000",
        pct_removed_hg38="77.78", pct_removed_t2t_phix="25.00",
        pct_removed_pangenome="33.33", pct_retained_overall="10.00",
        bracken_microbial_reads="250000", cleared_consensus_gate="true",
        bracken_genera="300", bracken_species="700",
        metaphlan_genera="90", metaphlan_species="180",
        shared_genera="75", shared_species="120",
        decontam_sample_type="Tumor", decontam_role="positive",
        decontam_final_status="retained",
    )
    row.update(over)
    return row


def render(rows, record=None):
    return CR.render(rows, record, Path("cmp.qc.summary.tsv"))


class SelfContainmentTests(unittest.TestCase):
    def test_no_network_dependency(self):
        """It must open from a filesystem, offline, years from now."""
        html = render([sample_row()])
        for forbidden in ("http://", "https://", "<script", "cdn."):
            self.assertNotIn(forbidden, html, f"page references {forbidden}")

    def test_has_dark_mode_and_viewport(self):
        html = render([sample_row()])
        self.assertIn("prefers-color-scheme:dark", html)
        self.assertIn('name="viewport"', html)


class ProvenanceTests(unittest.TestCase):
    def test_dirty_tree_is_stated_on_the_page(self):
        record = {"status": "completed", "environment": {"nextflow_version": "26.04.6"},
                  "run": {"started": "t0"},
                  "code": {"commit": "abc123456789", "dirty": True,
                           "modified_files": ["M main.nf", "M x.nf"]}}
        html = render([sample_row()], record)
        self.assertIn("Uncommitted changes", html)
        self.assertIn("2 modified path", html)

    def test_clean_tree_has_no_warning(self):
        record = {"status": "completed", "environment": {}, "run": {},
                  "code": {"commit": "abc", "dirty": False}}
        self.assertNotIn("Uncommitted changes", render([sample_row()], record))

    def test_absent_provenance_is_called_out(self):
        html = render([sample_row()], None)
        self.assertIn("No provenance record", html)


class ChartGeometryTests(unittest.TestCase):
    def geometry(self, html):
        vb = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', html)
        rects = re.findall(
            r'class="seg" x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="([\d.]+)"', html)
        return (float(vb.group(1)), float(vb.group(2))), [tuple(map(float, r)) for r in rects]

    def test_nothing_overflows_the_viewbox(self):
        rows = [sample_row(f"S{i}") for i in range(6)]
        (w, h), rects = self.geometry(render(rows))
        self.assertTrue(rects, "no segments drawn")
        for x, y, rw, rh in rects:
            self.assertLessEqual(x + rw, w + 0.01, "segment overflows right edge")
            self.assertLessEqual(y + rh, h + 0.01, "segment overflows bottom edge")
            self.assertGreaterEqual(rw, 0.0, "negative width")

    def test_row_count_matches_sample_count(self):
        rows = [sample_row(f"S{i}") for i in range(5)]
        _, rects = self.geometry(render(rows))
        self.assertEqual(len({r[1] for r in rects}), 5)

    def test_identity_is_never_colour_alone(self):
        """A legend for the four series, plus the same numbers in a table."""
        html = render([sample_row()])
        self.assertEqual(html.count('class="key"'), 4)
        self.assertIn("Filtering and depletion", html)
        self.assertIn("<title>", html)      # per-segment hover


class LibraryDenominatorTests(unittest.TestCase):
    """Audit C07: BAM/CRAM retention must not read as retention of the whole library."""

    def test_bam_cohort_shows_library_tiles(self):
        row = sample_row(library_primary_records="31", pct_retained_of_library="90.32")
        html = render([row])
        self.assertIn("Library records (BAM/CRAM)", html)
        self.assertIn("Median retention of library", html)
        self.assertIn("only the unmapped reads extracted", html)
        self.assertNotIn("before filtering, all samples", html)   # the old, misleading note

    def test_fastq_cohort_has_no_library_tiles(self):
        self.assertNotIn("Library records (BAM/CRAM)", render([sample_row()]))


class ValidationSectionTests(unittest.TestCase):
    def test_unscored_cells_are_counted(self):
        """A15: kept-but-unvalidated cells must be visible in the report."""
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "validation_mask_report.tsv"
            p.write_text("clade_name\tn_validated\tn_failed\tn_insufficient\tn_not_in_reference\tn_absent\t"
                         "n_unscored\treads_removed\nd__B|g__G|s__G_a\t2\t1\t0\t0\t3\t4\t50\n")
            html = CR.validation_section(p)
            self.assertIn("<b>4</b> unscored", html)


class DegradationTests(unittest.TestCase):
    def test_empty_cohort_does_not_crash(self):
        html = render([])
        self.assertIn("CMPipeline cohort report", html)

    def test_missing_depletion_explains_itself(self):
        row = sample_row()
        for key in ("pct_retained_overall", "pct_removed_hg38",
                    "pct_removed_t2t_phix", "pct_removed_pangenome"):
            row[key] = "NA"
        html = render([row])
        self.assertIn("No depletion counts", html)

    def test_sample_below_gate_shows_dashes_not_zeros(self):
        row = sample_row(cleared_consensus_gate="false", metaphlan_genera="NA",
                         metaphlan_species="NA", shared_genera="NA", shared_species="NA")
        html = render([row])
        self.assertIn("Cleared gate", html)
        self.assertIn("—", html)

    def test_decontam_section_absent_when_it_did_not_run(self):
        row = sample_row(decontam_role="NA", decontam_sample_type="NA",
                         decontam_final_status="NA")
        self.assertNotIn("Decontamination roles", render([row]))


class EndToEndTests(unittest.TestCase):
    def test_writes_a_file_from_a_real_tsv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            qc = tmp / "cmp.qc.summary.tsv"
            with qc.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=QC.COLUMNS, delimiter="\t",
                                        extrasaction="ignore")
                writer.writeheader()
                writer.writerows([sample_row("A"), sample_row("B")])
            prov = tmp / "provenance"; prov.mkdir()
            (prov / "s.json").write_text(json.dumps(
                {"status": "completed", "code": {"commit": "abc", "dirty": False},
                 "environment": {}, "run": {}}))
            out = tmp / "report.html"
            import sys
            argv = sys.argv
            sys.argv = ["build_cohort_report.py", "--qc-summary", str(qc),
                        "--output", str(out), "--provenance", str(prov)]
            try:
                CR.main()
            finally:
                sys.argv = argv
            self.assertTrue(out.is_file())
            text = out.read_text()
            self.assertIn("A", text)
            self.assertIn("B", text)


if __name__ == "__main__":
    unittest.main()
