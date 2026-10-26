#!/usr/bin/env python3
"""QC floor: a sample that cannot be assessed is `not_assessable`, never a microbial negative.

pksProfiler COHORT_CHECKS A3 (a proposal there; implemented here), port review #34. Three reasons:
  no_unmapped_records              the BAM/CRAM has zero FLAG 4 records
  unmapped_filtered_upstream       the input's @PG chain shows samtools view -F 4 (etc.) upstream
  below_min_reads_after_depletion  fewer host-depleted reads than --min_reads_after_depletion (> 0)

The extraction tests run Modules/extract_reads.nf's rendered script on synthetic BAMs (helpers from
tests/test_low_score_rescue.py); everything else uses hand-written synthetic QC artefacts.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_low_score_rescue as LSRT   # noqa: E402  (synthetic SAM records + the module renderer)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


QC = load("build_qc_summary_floor", ROOT / "scripts" / "build_qc_summary.py")
FP = load("finalize_provenance_floor", ROOT / "scripts" / "finalize_provenance.py")
CR = load("build_cohort_report_floor", ROOT / "scripts" / "build_cohort_report.py")
MAIN = (ROOT / "main.nf").read_text()
MODULE = (ROOT / "Modules" / "extract_reads.nf").read_text()

PG_BWA = "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem -t 8 ref.fa SECRET_1.fq SECRET_2.fq"
PG_FILTER = ("@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.9\t"
             "CL:samtools view -b -F 4 -o /SECRET/mapped_only.bam /SECRET/in.bam")
MAPPED_ONLY = [r for r in LSRT.RECORDS if r[2]]          # no FLAG 4 record at all


@unittest.skipUnless(LSRT.SAMTOOLS_BIN, "samtools >= 1.16 + htsfile + bgzip not found (set CMP_TEST_SAMTOOLS_BIN)")
class ExtractionRecordsTheFloor(unittest.TestCase):
    def run_module(self, pg, records=LSRT.RECORDS, mode="auto"):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        (tmp / "in.sam").write_text(LSRT.sam_text(pg, records))
        subprocess.run([str(LSRT.SAMTOOLS_BIN / "samtools"), "view", "-b", "-o", str(tmp / "in.bam"),
                        str(tmp / "in.sam")], check=True)
        work = tmp / "task"
        work.mkdir()
        script = LSRT.render_script(mode, str(tmp / "in.bam"))
        # Nextflow's env(INPUT_QC_FLOOR) output reads the variable after the script; do the same.
        (work / ".command.sh").write_text(script + '\necho "$INPUT_QC_FLOOR" > .floor\n')
        env = {"PATH": f"{LSRT.SAMTOOLS_BIN}:/usr/bin:/bin", "HOME": str(tmp), "LC_ALL": "C"}
        r = subprocess.run(["bash", ".command.sh"], cwd=work, env=env, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        counts = LSRT.read_tsv_row(work / "S.extract_counts.tsv")
        return counts, (work / ".floor").read_text().strip(), r.stderr

    def test_normal_input_passes(self):
        counts, floor, _ = self.run_module((PG_BWA,))
        self.assertEqual((counts["input_qc_floor"], floor), ("pass", "pass"))

    def test_upstream_unmapped_filter(self):
        counts, floor, err = self.run_module((PG_BWA, PG_FILTER), MAPPED_ONLY)
        self.assertEqual((counts["input_qc_floor"], floor), ("unmapped_filtered_upstream",) * 2)
        self.assertIn("NOT ASSESSABLE", err)
        self.assertNotIn("SECRET", err)

    def test_upstream_filter_wins_even_with_unmapped_records_present(self):
        counts, floor, _ = self.run_module((PG_BWA, PG_FILTER))
        self.assertEqual(floor, "unmapped_filtered_upstream")

    def test_no_flag4_records(self):
        for mode in ("off", "on"):              # the rescue adding low-AS reads does not lift it
            with self.subTest(mode=mode):
                counts, floor, _ = self.run_module((PG_BWA,), MAPPED_ONLY, mode=mode)
                self.assertEqual(counts["extracted_unmapped_records"], "0")
                self.assertEqual(floor, "no_unmapped_records")


def write_row(path: Path, header: list[str], values: list) -> None:
    path.write_text("\t".join(header) + "\n" + "\t".join(str(v) for v in values) + "\n")


def sample_inputs(d: Path, sample: str, unmapped: int, upstream: str = "none", floor: str | None = None,
                  depleted: int | None = 500, fastp: int = 1000, input_type: str = "bam"):
    """Synthetic per-sample QC artefacts for one BAM sample."""
    if floor is None:
        floor = ("unmapped_filtered_upstream" if upstream != "none" else
                 "no_unmapped_records" if unmapped == 0 else "pass")
    write_row(d / f"{sample}.extract_counts.tsv",
              ["sample", "library_primary_records", "extracted_unmapped_records", "low_as_reads_extracted",
               "low_as_no_as_records", "input_qc_floor"], [sample, 10000, unmapped, 0, "NA", floor])
    write_row(d / f"{sample}.input_provenance.tsv", ["sample", "input_ref_fingerprint", "input_unmapped_filter"],
              [sample, "a" * 32, upstream])
    if floor in ("unmapped_filtered_upstream", "no_unmapped_records"):
        return                                     # --qc_floor_action exclude: nothing downstream ran
    (d / f"{sample}.fastp.json").write_text(json.dumps({"summary": {
        "before_filtering": {"total_reads": unmapped}, "after_filtering": {"total_reads": fastp}}}))
    if depleted is not None:
        (d / f"{sample}.depletion_counts.tsv").write_text(
            "sample\tstage\treads\n" + "".join(f"{sample}\t{st}\t{n}\n" for st, n in (
                ("input", fastp), ("after_hg38", depleted), ("after_t2t_phix", depleted),
                ("host_depleted", depleted))))


class QcSummaryTerminalState(unittest.TestCase):
    STAGES = ["extraction", "read_qc", "host_depletion"]

    def build(self, make, min_after=0, types=None):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make(d)
            samples = sorted(p.name.split(".")[0] for p in d.glob("*.extract_counts.tsv"))
            (d / "sample_input_types.tsv").write_text("".join(f"{s}\t{(types or {}).get(s, 'bam')}\n"
                                                              for s in samples))
            (d / "expected_samples.txt").write_text("\n".join(samples) + "\n")
            return {r["sample"]: r for r in QC.build(d, 100000, self.STAGES, min_after)}

    def test_reasons_and_precedence(self):
        def make(d):
            sample_inputs(d, "ok", 1000)
            sample_inputs(d, "filtered", 0, upstream="samtools_view_exclude_unmapped")
            sample_inputs(d, "empty", 0)
            sample_inputs(d, "thin", 1000, depleted=10)
        rows = self.build(make)
        self.assertEqual(rows["ok"]["terminal_state"], "completed")
        self.assertEqual(rows["filtered"]["terminal_state"], "not_assessable_unmapped_filtered_upstream")
        self.assertEqual(rows["empty"]["terminal_state"], "not_assessable_no_unmapped_records")
        self.assertEqual(rows["thin"]["terminal_state"], "completed")          # floor off by default
        for s in ("filtered", "empty"):
            # Excluded at extraction: no read_qc output, but that is the floor, not a failure.
            self.assertNotIn("failed_at", rows[s]["terminal_state"])
            self.assertEqual(rows[s]["sample_status"], "not_assessable")
        self.assertEqual(rows["ok"]["qc_floor"], "pass")

    def test_depth_floor(self):
        def make(d):
            sample_inputs(d, "ok", 1000, depleted=500)
            sample_inputs(d, "thin", 1000, depleted=10)
            sample_inputs(d, "zero_fastp", 1000, fastp=0, depleted=0)
        rows = self.build(make, min_after=100)
        self.assertEqual(rows["ok"]["qc_floor"], "pass")
        self.assertEqual(rows["thin"]["terminal_state"], "not_assessable_below_min_reads_after_depletion")
        self.assertEqual(rows["thin"]["sample_status"], "not_assessable")
        self.assertEqual(rows["zero_fastp"]["terminal_state"], "excluded_zero_reads_after_fastp")

    def test_fastq_input_is_judged_on_depth_only(self):
        def make(d):
            (d / "fq.fastp.json").write_text(json.dumps({"summary": {
                "before_filtering": {"total_reads": 50}, "after_filtering": {"total_reads": 50}}}))
            (d / "fq.depletion_counts.tsv").write_text(
                "sample\tstage\treads\nfq\tinput\t50\nfq\tafter_hg38\t5\nfq\tafter_t2t_phix\t5\nfq\thost_depleted\t5\n")
            (d / "fq.extract_counts.tsv").unlink(missing_ok=True)
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make(d)
            (d / "sample_input_types.tsv").write_text("fq\tfastq\n")
            off = QC.build(d, 100000, ["read_qc", "host_depletion"], 0)[0]
            on = QC.build(d, 100000, ["read_qc", "host_depletion"], 10)[0]
        self.assertEqual(off["qc_floor"], "NA")          # nothing to judge a FASTQ on with the floor off
        self.assertEqual(on["terminal_state"], "not_assessable_below_min_reads_after_depletion")

    def test_recorded_and_recomputed_floor_must_agree(self):
        def make(d):
            sample_inputs(d, "liar", 1000, floor="no_unmapped_records")
        with self.assertRaisesRegex(ValueError, "input_qc_floor"):
            self.build(make)

    def test_failure_at_extraction_is_still_a_failure(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "sample_input_types.tsv").write_text("gone\tbam\n")
            (d / "expected_samples.txt").write_text("gone\n")
            row = QC.build(d, 100000, self.STAGES)[0]
        self.assertEqual(row["terminal_state"], "failed_at_extraction")
        self.assertEqual(row["qc_floor"], "NA")


class Reporting(unittest.TestCase):
    def qc_dir(self, d: Path, states):
        out = d / "QC"
        out.mkdir()
        with open(out / "cmp.qc.summary.tsv", "w") as fh:
            fh.write("sample\tterminal_state\n")
            for i, s in enumerate(states):
                fh.write(f"SECRET_{i}\t{s}\n")
        return out

    def test_run_report_note_counts_only(self):
        with tempfile.TemporaryDirectory() as d:
            out = self.qc_dir(Path(d), ["completed_profiled", "not_assessable_no_unmapped_records",
                                        "not_assessable_unmapped_filtered_upstream"])
            outcomes = FP.sample_outcomes({"qc_summary_dir": str(out)})
            self.assertEqual(outcomes["not_assessable"], 2)
            text = "\n".join(FP.render_sample_outcomes(outcomes))
            self.assertIn("2 sample(s) NOT ASSESSABLE", text)
            self.assertIn("NOT microbial negatives", text)
            self.assertIn("in no cohort table", text)
            flagged = FP.sample_outcomes({"qc_summary_dir": str(out), "qc_floor_action": "flag"})
            self.assertIn("WERE profiled", "\n".join(FP.render_sample_outcomes(flagged)))
            self.assertNotIn("SECRET", text)

    def test_cohort_report_tile(self):
        rows = [{"terminal_state": "completed_profiled"}, {"terminal_state": "not_assessable_no_unmapped_records"}]
        self.assertIn("Not assessable (QC floor)", CR.stat_tiles(rows))
        self.assertNotIn("Not assessable", CR.stat_tiles(rows[:1]))


class Wiring(unittest.TestCase):
    LIB = (ROOT / "lib" / "QcFloor.groovy").read_text()

    def test_params_and_checks(self):
        self.assertIn('params.qc_floor_action = "exclude"', MAIN)
        self.assertIn("params.min_reads_after_depletion = 0", MAIN)
        self.assertIn("def qcFloorProblem = QcFloor.problem(params)", MAIN)
        self.assertIn("if (qcFloorProblem) error qcFloorProblem", MAIN)
        self.assertIn("--qc_floor_action must be exclude or flag", self.LIB)
        self.assertIn("--min_reads_after_depletion must be a whole number >= 0", self.LIB)

    def test_extraction_gate(self):
        self.assertIn("tuple val(sampleID), env('INPUT_QC_FLOOR'), emit: floor", MODULE)
        self.assertIn("extractReads.out.reads.join(extractReads.out.floor)", MAIN)
        self.assertIn("keep: QcFloor.keepExtracted(floor, qcFloorAction)", MAIN)
        # The gated channel, not the raw one, feeds read QC.
        self.assertNotIn("extractReads.out.reads.set { EXTRACTED_READS }", MAIN)
        self.assertIn(".mix(EXTRACTED_READS.map", MAIN)

    def test_depth_gate_only_when_on(self):
        self.assertIn("if (minReadsAfterDepletion > 0) {", MAIN)
        self.assertIn("MAPPED_READS_MULTI.PAN.set { CLASSIFICATION_READS }", MAIN)
        self.assertIn("keep: QcFloor.keepDepleted(n, minReadsAfterDepletion, qcFloorAction)", MAIN)
        qc = (ROOT / "Modules" / "qc_summary.nf").read_text()
        self.assertIn('--min-reads-after-depletion "${params.min_reads_after_depletion}"', qc)

    @unittest.skipUnless(shutil.which("nextflow") and shutil.which("java"), "needs nextflow and java")
    def test_gates_execute(self):
        """The exact channel code from main.nf, with lib/QcFloor.groovy, on toy channels."""
        ext = MAIN.split("            extractReads.out.reads.join(extractReads.out.floor)", 1)[1]
        ext = "extractReads.out.reads.join(extractReads.out.floor)" + ext.split("            qcParts <<", 1)[0]
        dep = MAIN.split("            if (minReadsAfterDepletion > 0) {", 1)[1].split("            } else {", 1)[0]
        for action, expect_ext, expect_dep in (("exclude", "[ok]", "[deep]"), ("flag", "[filt, ok]", "[deep, thin]")):
            with self.subTest(action=action), tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                (tmp / "lib").mkdir()
                shutil.copy(ROOT / "lib" / "QcFloor.groovy", tmp / "lib")
                (tmp / "deep.tsv").write_text("sample\tstage\treads\ndeep\thost_depleted\t500\n")
                (tmp / "thin.tsv").write_text("sample\tstage\treads\nthin\tinput\t90\nthin\thost_depleted\t5\n")
                (tmp / "t.nf").write_text(
                    "workflow {\n"
                    f"    def qcFloorAction = '{action}'\n"
                    "    def minReadsAfterDepletion = 100L\n"
                    "    def log = [warn: { m -> println m }]\n"
                    "    def extractReads = [out: [reads: Channel.of(tuple('ok', 'r1'), tuple('filt', 'r2')),\n"
                    "                              floor: Channel.of(tuple('filt', 'unmapped_filtered_upstream'), tuple('ok', 'pass'))]]\n"
                    f"    {ext}\n"
                    "    EXTRACTED_READS.map { it[0] }.toSortedList().view { \"EXTRACTED ${it}\" }\n"
                    "    def MAPPED_READS_MULTI = [PAN: Channel.of(tuple('deep', 'd'), tuple('thin', 't'))]\n"
                    "    def mapReads = [out: [counts: Channel.of(file(\"${launchDir}/deep.tsv\"), file(\"${launchDir}/thin.tsv\"))]]\n"
                    f"    {dep}\n"
                    "    CLASSIFICATION_READS.map { it[0] }.toSortedList().view { \"CLASSIFIED ${it}\" }\n"
                    "}\n")
                env = {k: v for k, v in os.environ.items() if k != "LD_LIBRARY_PATH"}
                r = subprocess.run(["nextflow", "-q", "run", "t.nf"], cwd=tmp, env=env,
                                   capture_output=True, text=True, timeout=600)
                self.assertEqual(r.returncode, 0, r.stdout[-3000:] + r.stderr[-3000:])
                self.assertIn(f"EXTRACTED {expect_ext}", r.stdout)
                self.assertIn(f"CLASSIFIED {expect_dep}", r.stdout)
                if action == "exclude":
                    self.assertIn("Sample filt NOT ASSESSABLE (unmapped_filtered_upstream)", r.stdout + r.stderr)
                    self.assertIn("Sample thin NOT ASSESSABLE (5 reads after host depletion)", r.stdout + r.stderr)

    def test_groovy_escapes_in_the_floor_awk(self):
        self.assertIn("awk -F '\\\\t' 'NR == 1", MODULE)

    def test_documented(self):
        outputs = (ROOT / "docs" / "outputs.md").read_text()
        params = (ROOT / "docs" / "parameters.md").read_text()
        for word in ("not_assessable_no_unmapped_records", "not_assessable_unmapped_filtered_upstream",
                     "not_assessable_below_min_reads_after_depletion", "qc_floor"):
            self.assertIn(word, outputs)
        for p in ("--qc_floor_action", "--min_reads_after_depletion"):
            self.assertIn(p, params)
        self.assertIn("not_assessable", (ROOT / "docs" / "CHANGELOG.md").read_text())


if __name__ == "__main__":
    unittest.main()
