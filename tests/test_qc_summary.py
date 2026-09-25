#!/usr/bin/env python3
"""Checks for scripts/build_qc_summary.py, on synthetic fixtures only.

Covers the things that would silently produce a wrong table: the consensus gate must
use the same definition as Bracken.nf, retention percentages must chain stage to stage,
a stage that did not run must yield NA rather than a zero, and the two profilers'
taxon names must be normalised to the same spelling before they are intersected.
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_qc_summary", ROOT / "scripts" / "build_qc_summary.py"
)
QC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QC)


def write_fastp(directory, sample, raw=1000, kept=900, adapter=True):
    payload = {
        "summary": {
            "before_filtering": {"total_reads": raw, "total_bases": raw * 100, "q30_rate": 0.95},
            "after_filtering": {"total_reads": kept, "total_bases": kept * 100, "q30_rate": 0.97},
        },
        "filtering_result": {"low_quality_reads": 60, "too_short_reads": 40},
        "duplication": {"rate": 0.11},
    }
    if adapter:
        payload["adapter_cutting"] = {"adapter_trimmed_reads": 25}
    (directory / f"{sample}.fastp.json").write_text(json.dumps(payload))


def write_depletion(directory, sample, values=(900, 500, 300, 200)):
    rows = zip(("input", "after_hg38", "after_t2t_phix", "host_depleted"), values)
    lines = ["sample\tstage\treads"] + [f"{sample}\t{stage}\t{n}" for stage, n in rows]
    (directory / f"{sample}.depletion_counts.tsv").write_text("\n".join(lines) + "\n")


def write_krakenuniq(directory, sample, root=250000, unclassified=50000):
    """A real KrakenUniq report: 9 columns, and BOTH rows carry rank "no rank".

    The earlier fixture here invented a 6-column layout with rank "U"/"R", which is why
    a parser that could never work against real output still passed its test. Real
    columns: % reads taxReads kmers dup cov taxID rank taxName
    """
    header = "# KrakenUniq v1.0.4 DB:/fake"
    lines = [
        header,
        f"5.0\t{unclassified}\t{unclassified}\t0\t0\t0\t0\tno rank\tunclassified",
        f"95.0\t{root}\t100\t5000\t1.1\t0.01\t1\tno rank\troot",
    ]
    (directory / f"{sample}.krakenuniq.report.txt").write_text("\n".join(lines) + "\n")


def write_bracken(directory, sample, level, names):
    header = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    lines = [header] + [f"{n}\t123\t{level}\t10\t5\t15\t0.1" for n in names]
    (directory / f"{sample}.bracken.{level}.report.txt").write_text("\n".join(lines) + "\n")


def write_metaphlan(directory, sample, clades):
    lines = ["#mpa_vJun23", "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional"]
    lines += [f"{c}\t1|2\t1.5\t" for c in clades]
    (directory / f"{sample}.profiled_metagenome.txt").write_text("\n".join(lines) + "\n")


def build_rows(directory, min_reads=100000):
    return {r["sample"]: r for r in QC.build(Path(directory), min_reads)}


class FastpTests(unittest.TestCase):
    def test_reads_and_pass_rate(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fastp(Path(tmp), "S1", raw=1000, kept=900)
            row = build_rows(tmp)["S1"]
            self.assertEqual(row["reads_raw"], 1000)
            self.assertEqual(row["reads_after_fastp"], 900)
            self.assertEqual(row["pct_pass_fastp"], "90.00")
            self.assertEqual(row["duplication_rate"], 0.11)

    def test_adapter_scan_recorded_both_ways(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fastp(Path(tmp), "WITH", adapter=True)
            write_fastp(Path(tmp), "WITHOUT", adapter=False)
            rows = build_rows(tmp)
            self.assertEqual(rows["WITH"]["adapter_scan_applied"], "true")
            self.assertEqual(rows["WITH"]["adapter_trimmed_reads"], 25)
            self.assertEqual(rows["WITHOUT"]["adapter_scan_applied"], "false")
            self.assertEqual(rows["WITHOUT"]["adapter_trimmed_reads"], "NA")


class DepletionTests(unittest.TestCase):
    def test_stage_percentages_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fastp(Path(tmp), "S1", raw=1000, kept=900)
            write_depletion(Path(tmp), "S1", values=(900, 450, 300, 150))
            row = build_rows(tmp)["S1"]
            self.assertEqual(row["pct_removed_hg38"], "50.00")        # 900 -> 450
            self.assertEqual(row["pct_removed_t2t_phix"], "33.33")    # 450 -> 300
            self.assertEqual(row["pct_removed_pangenome"], "50.00")   # 300 -> 150
            self.assertEqual(row["pct_retained_overall"], "15.00")    # 150 of 1000 raw

    def test_absent_depletion_is_na_not_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fastp(Path(tmp), "S1")
            row = build_rows(tmp)["S1"]
            for column in ("reads_host_depleted", "pct_removed_hg38", "pct_retained_overall"):
                self.assertEqual(row[column], "NA", column)


def write_extract_counts(directory, sample, library, unmapped):
    (directory / f"{sample}.extract_counts.tsv").write_text(
        "sample\tlibrary_primary_records\textracted_unmapped_records\n"
        f"{sample}\t{library}\t{unmapped}\n")


class LibraryDenominatorTests(unittest.TestCase):
    """Audit C07 / pksProfiler F09: for BAM/CRAM, reads_raw is the extracted unmapped subset."""

    def test_bam_sample_reports_retention_of_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_extract_counts(Path(tmp), "B", library=31, unmapped=28)
            write_fastp(Path(tmp), "B", raw=28, kept=28)
            write_depletion(Path(tmp), "B", values=(28, 28, 28, 28))
            row = build_rows(tmp)["B"]
            self.assertEqual(row["pct_retained_overall"], "100.00")      # of the extracted subset
            self.assertEqual(row["library_primary_records"], "31")
            self.assertEqual(row["pct_library_unmapped"], "90.32")       # 28 / 31
            self.assertEqual(row["pct_retained_of_library"], "90.32")    # 28 / 31, not 100

    def test_fastq_sample_has_na_library_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fastp(Path(tmp), "F")
            write_depletion(Path(tmp), "F")
            row = build_rows(tmp)["F"]
            for column in ("library_primary_records", "pct_library_unmapped", "pct_retained_of_library"):
                self.assertEqual(row[column], "NA", column)


class SampleStatusTests(unittest.TestCase):
    """Audit C13/C02: a no-call must not look like 'ran and found nothing'."""

    def test_status_per_outcome(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_fastp(d, "EMPTY", raw=100, kept=0)                       # fastp emptied it
            write_fastp(d, "THIN"); write_depletion(d, "THIN")
            (d / "THIN.bracken.G.report.txt").write_text(                 # header-only no-call
                "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n")
            write_fastp(d, "OK"); write_depletion(d, "OK"); write_bracken(d, "OK", "G", ["Gen"])
            rows = build_rows(tmp)
            self.assertEqual(rows["EMPTY"]["sample_status"], "dropped_zero_after_fastp")
            self.assertEqual(rows["EMPTY"]["reads_host_depleted"], 0)     # 0, not NA
            self.assertEqual(rows["THIN"]["bracken_status"], "no_call")
            self.assertEqual(rows["THIN"]["sample_status"], "below_bracken_threshold")
            self.assertEqual(rows["OK"]["sample_status"], "profiled")

    def test_sample_with_no_outputs_is_named(self):
        """C10: a sample whose tasks failed or were ignored must appear, not vanish."""
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_fastp(d, "OK")
            (d / "expected_samples.txt").write_text("OK\nGONE\n")
            rows = build_rows(tmp)
            self.assertEqual(rows["GONE"]["sample_status"], "no_outputs")
            self.assertNotEqual(rows["OK"]["sample_status"], "no_outputs")


class GateTests(unittest.TestCase):
    """The gate silently returned 0 for every sample from 2026-09-04 to 2026-09-23.

    Three compounding faults: the rank was read from column 4 (which is `kmers`), the
    rank of both rows is the literal string "no rank" so no rank comparison could match,
    and root/unclassified are siblings so subtracting one from the other was wrong
    arithmetic. Root's own clade count is the classified total.
    """

    def test_gate_is_root_clade_count(self):
        """Must match Bracken.nf exactly, or the table contradicts the pipeline."""
        with tempfile.TemporaryDirectory() as tmp:
            write_krakenuniq(Path(tmp), "S1", root=250000, unclassified=50000)
            row = build_rows(tmp, min_reads=100000)["S1"]
            self.assertEqual(row["bracken_microbial_reads"], 250000)
            self.assertEqual(row["cleared_consensus_gate"], "true")

    def test_unclassified_is_not_subtracted(self):
        """root and unclassified are siblings, not parent and child."""
        with tempfile.TemporaryDirectory() as tmp:
            write_krakenuniq(Path(tmp), "S1", root=250000, unclassified=900000)
            row = build_rows(tmp, min_reads=100000)["S1"]
            self.assertEqual(row["bracken_microbial_reads"], 250000)

    def test_below_threshold_is_false(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_krakenuniq(Path(tmp), "S1", root=42232, unclassified=85900)
            row = build_rows(tmp, min_reads=100000)["S1"]
            self.assertEqual(row["bracken_microbial_reads"], 42232)
            self.assertEqual(row["cleared_consensus_gate"], "false")

    def test_comment_lines_are_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_krakenuniq(Path(tmp), "S1", root=150000, unclassified=10)
            row = build_rows(tmp, min_reads=100000)["S1"]
            self.assertEqual(row["bracken_microbial_reads"], 150000)


class TaxonCountTests(unittest.TestCase):
    def test_counts_and_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_bracken(d, "S1", "G", ["Escherichia", "Bacteroides", "Fusobacterium"])
            write_bracken(d, "S1", "S", ["Escherichia coli", "Bacteroides fragilis"])
            write_metaphlan(d, "S1", [
                "k__Bacteria|g__Escherichia",
                "k__Bacteria|g__Bacteroides",
                "k__Bacteria|g__Akkermansia",
                "k__Bacteria|g__Escherichia|s__Escherichia_coli",
            ])
            row = build_rows(tmp)["S1"]
            self.assertEqual(row["bracken_genera"], 3)
            self.assertEqual(row["bracken_species"], 2)
            self.assertEqual(row["metaphlan_genera"], 3)
            self.assertEqual(row["metaphlan_species"], 1)
            # Escherichia and Bacteroides are shared; Fusobacterium and Akkermansia are not
            self.assertEqual(row["shared_genera"], 2)
            # s__Escherichia_coli must normalise to "Escherichia coli" to match Bracken
            self.assertEqual(row["shared_species"], 1)

    def test_metaphlan_absent_leaves_na(self):
        """A sample below the gate has no MetaPhlAn profile; that must read NA."""
        with tempfile.TemporaryDirectory() as tmp:
            write_bracken(Path(tmp), "S1", "G", ["Escherichia"])
            row = build_rows(tmp)["S1"]
            self.assertEqual(row["bracken_genera"], 1)
            self.assertEqual(row["metaphlan_genera"], "NA")
            self.assertEqual(row["shared_genera"], "NA")


class DecontamRoleTests(unittest.TestCase):
    def test_roles_are_joined(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_fastp(d, "S1")
            (d / "run.sample_validation.tsv").write_text(
                "sample_id\tin_otu\tsample_type\trole\tfinal_status\n"
                "S1\tTRUE\tTumor\tpositive\tretained\n"
            )
            row = build_rows(tmp)["S1"]
            self.assertEqual(row["decontam_role"], "positive")
            self.assertEqual(row["decontam_sample_type"], "Tumor")
            self.assertEqual(row["decontam_final_status"], "retained")


def write_sample_validation(directory, records):
    lines = ["sample_id\tin_otu\tin_metadata\tsample_type\trole\tfinal_status\treason"]
    lines += ["\t".join(r) for r in records]
    (directory / "decontam_run.sample_validation.tsv").write_text("\n".join(lines) + "\n")


def write_batch_ledger(directory, rank, records, sub="qc_1"):
    folder = directory / sub
    folder.mkdir(exist_ok=True)
    lines = ["sample\tstatus\treason"] + ["\t".join(r) for r in records]
    (folder / f"{rank}.batch_sample_status.tsv").write_text("\n".join(lines) + "\n")


class DownstreamStartTests(unittest.TestCase):
    """A run that starts at decontamination or batch correction still gets one row per sample."""

    def test_decontam_start_rows_come_from_sample_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_sample_validation(Path(tmp), [
                ("T1", "TRUE", "TRUE", "Primary Tumor", "positive", "included", "included_positive"),
                ("B1", "TRUE", "TRUE", "Blood Derived Normal", "control", "included", "included_control"),
                ("M1", "FALSE", "TRUE", "Primary Tumor", "NA", "excluded", "metadata_only"),
            ])
            rows = build_rows(tmp)
            self.assertEqual(sorted(rows), ["B1", "T1"])   # metadata-only is not an analysed sample
            self.assertEqual(rows["B1"]["decontam_role"], "control")
            self.assertEqual(rows["T1"]["reads_raw"], "NA")

    def test_metadata_only_row_still_joins_an_existing_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fastp(Path(tmp), "S1")
            write_sample_validation(Path(tmp), [
                ("S1", "FALSE", "TRUE", "Tumor", "NA", "excluded", "metadata_only")])
            self.assertEqual(build_rows(tmp)["S1"]["decontam_final_status"], "excluded")

    def test_batch_status_per_rank(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_batch_ledger(Path(tmp), "species", [
                ("T1", "included", "eligible"),
                ("B1", "excluded", "outside requested tumor subset"),
                ("X9", "unused", "metadata row absent from count table")])
            write_batch_ledger(Path(tmp), "genus", [
                ("T1", "included", "eligible"),
                ("B1", "excluded", "outside requested tumor subset")], sub="qc_2")
            rows = build_rows(tmp)
            self.assertEqual(sorted(rows), ["B1", "T1"])   # unused metadata row is not a sample
            self.assertEqual(rows["T1"]["batch_species_status"], "included")
            self.assertEqual(rows["T1"]["batch_genus_status"], "included")
            self.assertEqual(rows["B1"]["batch_genus_status"],
                             "excluded: outside requested tumor subset")

    def test_single_rank_leaves_the_other_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_batch_ledger(Path(tmp), "genus", [("T1", "included", "eligible")])
            row = build_rows(tmp)["T1"]
            self.assertEqual(row["batch_genus_status"], "included")
            self.assertEqual(row["batch_species_status"], "NA")

    def test_unknown_rank_ledger_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_batch_ledger(Path(tmp), "phylum", [("T1", "included", "eligible")])
            with self.assertRaises(ValueError):
                build_rows(tmp)


class StageConsistencyTests(unittest.TestCase):
    """Strict inputs, and a batch cohort that must be explained by decontamination's output."""

    DECONTAM = [
        ("T1", "TRUE", "TRUE", "Primary Tumor", "positive", "included", "included_positive"),
        ("T2", "TRUE", "TRUE", "Primary Tumor", "positive", "included", "included_positive"),
        ("B1", "TRUE", "TRUE", "Blood Derived Normal", "control", "included", "included_control"),
        ("M1", "FALSE", "TRUE", "Primary Tumor", "NA", "excluded", "metadata_only"),
    ]
    BATCH = [("T1", "included", "eligible"), ("T2", "excluded", "all-zero after decontamination"),
             ("B1", "excluded", "outside requested tumor subset")]

    def build(self, decontam=DECONTAM, species=BATCH, genus=BATCH, header=None):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = Path(tmp.name)
        if decontam is not None:
            write_sample_validation(d, decontam)
            if header is not None:
                path = d / "decontam_run.sample_validation.tsv"
                lines = path.read_text().splitlines()
                path.write_text("\n".join([header] + lines[1:]) + "\n")
        if species is not None:
            write_batch_ledger(d, "species", species)
        if genus is not None:
            write_batch_ledger(d, "genus", genus, sub="qc_2")
        return build_rows(d)

    def test_explained_exclusions_pass(self):
        rows = self.build()
        self.assertEqual(sorted(rows), ["B1", "T1", "T2"])
        self.assertEqual(rows["T2"]["batch_species_status"], "excluded: all-zero after decontamination")

    def test_genus_only_all_zero_exclusion_passes(self):
        genus = [("T1", "excluded", "all-zero after decontamination")] + self.BATCH[1:]
        self.assertEqual(self.build(genus=genus)["T1"]["batch_genus_status"],
                         "excluded: all-zero after decontamination")

    def test_missing_in_otu_column_is_an_error(self):
        with self.assertRaisesRegex(ValueError, "in_otu"):
            self.build(header="sample_id\tin_metadata\tsample_type\trole\tfinal_status\treason")

    def test_invalid_in_otu_value_is_an_error(self):
        for value in ("", "yes", "true", "NA"):
            with self.subTest(value=value):
                bad = [("T1", value) + self.DECONTAM[0][2:]] + self.DECONTAM[1:]
                with self.assertRaisesRegex(ValueError, "in_otu must be TRUE or FALSE"):
                    self.build(decontam=bad)

    def test_batch_sample_not_from_decontamination_is_an_error(self):
        extra = self.BATCH + [("X1", "included", "eligible")]
        with self.assertRaisesRegex(ValueError, "not from decontamination \\['X1'\\]"):
            self.build(species=extra, genus=extra)

    def test_decontaminated_sample_missing_from_batch_is_an_error(self):
        short = self.BATCH[:2]
        with self.assertRaisesRegex(ValueError, "missing from batch ledger \\['B1'\\]"):
            self.build(species=short, genus=short)

    def test_metadata_only_sample_in_batch_is_an_error(self):
        leaked = self.BATCH + [("M1", "excluded", "missing covariate: age_diag")]
        with self.assertRaisesRegex(ValueError, "not from decontamination \\['M1'\\]"):
            self.build(species=leaked, genus=leaked)

    def test_ranks_with_different_samples_are_an_error(self):
        with self.assertRaisesRegex(ValueError, "different samples per rank"):
            self.build(decontam=None, genus=self.BATCH[:2])

    def test_invalid_batch_ledger_rows_are_errors(self):
        for rows in ([("T1", "maybe", "eligible")], [("T1", "excluded", "")],
                     [("T1", "included", "eligible"), ("T1", "included", "eligible")],
                     [("T1", "unused", "metadata row absent from count table"), ("T1", "included", "eligible")]):
            with self.subTest(rows=rows):
                with self.assertRaises(ValueError):
                    self.build(decontam=None, species=rows, genus=None)

    def test_duplicate_decontam_sample_is_an_error(self):
        with self.assertRaisesRegex(ValueError, "duplicate sample T1"):
            self.build(decontam=self.DECONTAM + [self.DECONTAM[0]])


class OutputShapeTests(unittest.TestCase):
    def test_every_declared_column_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_fastp(Path(tmp), "S1")
            out = Path(tmp) / "qc.tsv"
            import sys
            argv = sys.argv
            sys.argv = ["build_qc_summary.py", "--inputs", tmp, "--output", str(out)]
            try:
                QC.main()
            finally:
                sys.argv = argv
            with out.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames, QC.COLUMNS)
                self.assertEqual(len(list(reader)), 1)


if __name__ == "__main__":
    unittest.main()
