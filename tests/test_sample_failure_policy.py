#!/usr/bin/env python3
"""Per-sample failure policy and the per-sample terminal state (audit X05/X06; pksProfiler dc5d703
and e126439).

X06: "Seven processes carry per_sample today. A new per-sample process added without the label
would silently get cohort-fatal finish." The per-sample processes are DISCOVERED from the source --
every process whose first input is `tuple val(sampleID), ...` -- rather than listed by hand
(pksProfiler's hand-kept map missed five module files, dc5d703), and each must carry
label 'per_sample', which conf/base.config gives --sample_failure_strategy.

X05: under --sample_failure_strategy ignore, (a) a sample whose Bracken task failed still had
host-depleted reads, so it was batched into alignForValidation, the scorer exited 1 on a sample
missing from the species table, and the run ended; (b) MetaPhlAn and HUMAnN failures were silent in
QC. Now validation batches only samples that are columns of the species table, and the QC summary
gives every sample exactly one terminal_state (completed_* / failed_at_<stage> / excluded_<reason>)
plus metaphlan_status and humann_status. The QC tests below build synthetic per-sample artefacts
(counts and taxon names only) and run scripts/build_qc_summary.py on them.
"""

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

REPO = Path(__file__).resolve().parents[1]
MODULES = sorted((REPO / "Modules").glob("*.nf"))
sys.path.insert(0, str(REPO / "scripts"))
SPEC = importlib.util.spec_from_file_location("build_qc_summary", REPO / "scripts" / "build_qc_summary.py")
QC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QC)

PROCESS_RE = re.compile(r"^process\s+(\w+)\s*\{", re.M)
FIRST_INPUT_RE = re.compile(r"\binput:\s*(?://[^\n]*\n|\s)*(.+)")
PER_SAMPLE_INPUT = re.compile(r"^(tuple\s+)?val\(\s*sampleID\s*\)")

# Processes whose first input is per-sample but which deliberately stay cohort-strict. Each needs a
# reason; empty on purpose (pksProfiler's list is empty too).
NOT_PER_SAMPLE_OK: dict = {}
# The per-sample processes as of 2026-10-26, so a discovery regression (finding nothing) fails.
KNOWN_PER_SAMPLE = {"extractReads", "FASTQC", "filterReads", "mapReads", "Bracken", "metaphlan4",
                    "humann3", "lowAsSurvival"}


def processes(text):
    starts = list(PROCESS_RE.finditer(text))
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        yield match.group(1), text[match.start():end]


def first_input(body):
    match = FIRST_INPUT_RE.search(body)
    return match.group(1).strip() if match else ""


def per_sample_processes(sources=None):
    """{name: has the per_sample label} for every process whose first input is per-sample."""
    found = {}
    for module in MODULES if sources is None else []:
        for name, body in processes(module.read_text()):
            if PER_SAMPLE_INPUT.match(first_input(body)):
                found[name] = "label 'per_sample'" in body
    for name, body in (sources or {}).items():
        if PER_SAMPLE_INPUT.match(first_input(body)):
            found[name] = "label 'per_sample'" in body
    return found


class PerSampleLabel(unittest.TestCase):

    def test_discovery_finds_the_known_processes(self):
        self.assertTrue(KNOWN_PER_SAMPLE <= set(per_sample_processes()),
                        sorted(KNOWN_PER_SAMPLE - set(per_sample_processes())))

    def test_every_per_sample_process_has_the_label(self):
        missing = sorted(n for n, labelled in per_sample_processes().items()
                         if not labelled and n not in NOT_PER_SAMPLE_OK)
        self.assertEqual(missing, [], "per-sample processes without label 'per_sample' "
                                      "(they would stop the cohort under --sample_failure_strategy ignore)")

    def test_only_per_sample_processes_have_the_label(self):
        """A cohort reducer must stay strict: the label on one would let it fail silently."""
        discovered = per_sample_processes()
        for module in MODULES:
            for name, body in processes(module.read_text()):
                if "label 'per_sample'" in body:
                    self.assertIn(name, discovered, f"{name} is labelled per_sample but its first input "
                                                    f"is not per-sample: {first_input(body)[:60]}")

    def test_removing_a_label_is_caught(self):
        """The check is not vacuous: the same discovery on a module with the label removed fails."""
        text = (REPO / "Modules" / "Bracken.nf").read_text()
        stripped = {name: body.replace("label 'per_sample'", "") for name, body in processes(text)}
        self.assertEqual(per_sample_processes(stripped), {"Bracken": False})

    def test_validation_batches_are_not_per_sample(self):
        """alignForValidation takes a batch (tuple val(sample_ids), ...): cohort-strict by design."""
        self.assertNotIn("alignForValidation", per_sample_processes())

    def test_allow_list_is_not_stale(self):
        self.assertEqual(sorted(set(NOT_PER_SAMPLE_OK) - set(per_sample_processes())), [])

    def test_label_wires_the_param(self):
        base = (REPO / "conf" / "base.config").read_text()
        block = base.split("withLabel: per_sample", 1)[1].split("}", 2)
        self.assertIn("params.sample_failure_strategy", "".join(block[:2]))
        main = (REPO / "main.nf").read_text()
        self.assertIn("params.sample_failure_strategy = 'finish'", main)
        self.assertIn("['finish', 'ignore']", main)


class DownstreamWiring(unittest.TestCase):
    """X05: downstream steps handle a sample that failed upstream explicitly."""

    @classmethod
    def setUpClass(cls):
        cls.main = (REPO / "main.nf").read_text()

    def test_validation_batches_only_samples_in_the_species_table(self):
        validation = self.main.split("STEP 3b: ALIGNMENT VALIDATION", 1)[1].split("STEP 4", 1)[0]
        self.assertIn("VALIDATION_TABLE_SAMPLES", validation)
        batches = validation.split("VALIDATION_BATCHES", 1)[0]
        self.assertIn("in_table", batches)
        self.assertRegex(batches, r"absent\.subscribe")   # logged, never silent
        self.assertNotRegex(batches, r"CLASSIFICATION_READS\s*\n\s*\.toSortedList")

    @unittest.skipUnless(shutil.which("nextflow") and shutil.which("java"), "needs nextflow and java")
    def test_validation_batch_filter_executes(self):
        """The exact channel code from main.nf, run on a toy table: B (not a column) is dropped and
        logged, A and C are batched in sorted order."""
        block = self.main.split("            VALIDATION_TABLE.first()", 1)[1]
        block = "VALIDATION_TABLE.first()" + block.split(".flatMap { all ->", 1)[0]
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "table.tsv").write_text("clade_name\tA\tC\nx\t1\t2\n")
            (Path(tmp) / "t.nf").write_text(
                "workflow {\n"
                "    Channel.of(file(\"${launchDir}/table.tsv\")).set { VALIDATION_TABLE }\n"
                "    Channel.of(tuple('C', 'rc'), tuple('B', 'rb'), tuple('A', 'ra')).set { CLASSIFICATION_READS }\n"
                "    def log = [warn: { m -> println m }]\n"
                f"    {block}\n"
                "        .view { \"BATCHED ${it}\" }\n"
                "}\n")
            env = {k: v for k, v in os.environ.items() if k != "LD_LIBRARY_PATH"}
            result = subprocess.run(["nextflow", "-q", "run", "t.nf"], cwd=tmp, env=env,
                                    capture_output=True, text=True, timeout=600)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("BATCHED [[A, ra], [C, rc]]", result.stdout)
            self.assertIn("sample B is not in the species table", result.stdout + result.stderr)

    def test_qc_summary_gets_stages_input_types_and_humann_completion(self):
        self.assertIn("sample_input_types.tsv", self.main)
        self.assertIn("humann3_completed.txt", self.main)
        self.assertRegex(self.main, r"qcSummary\(qcChannel\.collect\(\), qcStages,")
        module = (REPO / "Modules" / "qc_summary.nf").read_text()
        self.assertIn("val stages", module)
        self.assertIn('--stages "${stages}"', module)

    def test_stage_names_match_the_script(self):
        names = set(re.findall(r"'(extraction|read_qc|host_depletion|classification|metaphlan|humann|"
                               r"low_as_survival)'", self.main))
        self.assertEqual(names, set(QC.STAGE_ORDER))

    def test_default_stays_finish_and_is_documented(self):
        docs = "\n".join(p.read_text() for p in [REPO / "README.md", *sorted((REPO / "docs").rglob("*.md"))])
        self.assertIn("terminal_state", docs)
        self.assertIn("failed_at_", docs)
        self.assertRegex(docs, r"--sample_failure_strategy ignore")


# ------------------------------------------------------------------------- QC terminal state

def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def fastp(d, s, reads):
    write(d / f"{s}.fastp.json", json.dumps({"summary": {"before_filtering": {"total_reads": reads + 10},
                                                         "after_filtering": {"total_reads": reads}}}))


def depletion(d, s, reads):
    write(d / f"{s}.depletion_counts.tsv", "sample\tstage\treads\n" + "".join(
        f"{s}\t{st}\t{reads}\n" for st in ("input", "after_hg38", "after_t2t_phix", "host_depleted")))


def kraken(d, s, root):
    write(d / f"{s}.krakenuniq.report.txt",
          "#comment\n%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n"
          f"10\t5\t5\t0\t0\t0\t0\tno rank\tunclassified\n90\t{root}\t0\t0\t0\t0\t1\tno rank\troot\n")


def bracken(d, s, called=True):
    body = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"
    if called:
        body += "Escherichia\t561\tG\t50\t0\t50\t1.0\n"
    write(d / f"{s}.bracken.G.report.txt", body)
    write(d / f"{s}.bracken.S.report.txt", body.replace("Escherichia\t561\tG", "Escherichia coli\t562\tS"))


def metaphlan(d, s):
    write(d / f"{s}.profiled_metagenome.txt", "#mpa\nk__Bacteria|p__P|c__C|o__O|f__F|g__Escherichia\t100\n")


def extract(d, s):
    write(d / f"{s}.extract_counts.tsv",
          "sample\tlibrary_primary_records\textracted_unmapped_records\tlow_as_reads_extracted\tlow_as_no_as_records\n"
          f"{s}\t1000\t100\t0\tNA\n")


ALL = "extraction,read_qc,host_depletion,classification,metaphlan,humann,low_as_survival"


class TerminalState(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="cmp_terminal_"))
        self.d = self.tmp / "qc_1"

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def complete(self, s, root=200_000, mpa=True, humann=True, called=True):
        fastp(self.d, s, 1000)
        depletion(self.d, s, 900)
        kraken(self.d, s, root)
        bracken(self.d, s, called)
        if mpa:
            metaphlan(self.d, s)
        if humann:
            self.humann.append(s)

    def run_summary(self, samples, stages=ALL, types=None):
        write(self.d / "expected_samples.txt", "\n".join(samples) + "\n")
        write(self.d / "sample_input_types.tsv",
              "".join(f"{s}\t{(types or {}).get(s, 'fastq')}\n" for s in samples))
        if self.humann:
            write(self.d / "humann3_completed.txt", "\n".join(self.humann) + "\n")
        out = self.tmp / "cmp.qc.summary.tsv"
        args = [sys.executable, str(REPO / "scripts" / "build_qc_summary.py"), "--inputs", str(self.tmp),
                "--output", str(out), "--min-reads", "100000"]
        if stages is not None:
            args += ["--stages", stages]
        result = subprocess.run(args, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with out.open() as handle:
            return {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}, result.stdout

    humann: list

    def test_one_state_per_outcome(self):
        """Six samples, six outcomes -- the X05 forced-failure scenario, synthetic."""
        self.humann = []
        self.complete("OK")                                   # everything ran
        self.complete("THIN", root=10, mpa=False, humann=False, called=False)   # below gate, no call
        self.complete("NOMPA", mpa=False, humann=False)       # cleared the gate, MetaPhlAn failed
        self.complete("NOHUM", humann=False)                  # MetaPhlAn ran, HUMAnN failed
        fastp(self.d, "NOBRK", 1000); depletion(self.d, "NOBRK", 900)    # Bracken failed
        fastp(self.d, "EMPTY", 0)                             # fastp emptied it: dropped on purpose
        rows, stdout = self.run_summary(["OK", "THIN", "NOMPA", "NOHUM", "NOBRK", "EMPTY", "GONE"])
        state = {s: r["terminal_state"] for s, r in rows.items()}
        self.assertEqual(state, {
            "OK": "completed_profiled", "THIN": "completed_no_call", "NOMPA": "failed_at_metaphlan",
            "NOHUM": "failed_at_humann", "NOBRK": "failed_at_classification",
            "EMPTY": "excluded_zero_reads_after_fastp", "GONE": "failed_at_read_qc"})
        self.assertEqual({s: rows[s]["metaphlan_status"] for s in ("OK", "THIN", "NOMPA", "NOHUM", "NOBRK")},
                         {"OK": "ran", "THIN": "below_gate", "NOMPA": "failed", "NOHUM": "ran", "NOBRK": "NA"})
        self.assertEqual({s: rows[s]["humann_status"] for s in ("OK", "THIN", "NOMPA", "NOHUM")},
                         {"OK": "ran", "THIN": "below_gate", "NOMPA": "metaphlan_failed", "NOHUM": "failed"})
        self.assertIn("WARNING: 4 sample(s) failed", stdout)
        self.assertIn("failed_at_metaphlan: 1", stdout)

    def test_extraction_is_expected_of_alignments_only(self):
        self.humann = []
        self.complete("BAMOK"); extract(self.d, "BAMOK")
        self.complete("FQOK")
        rows, _ = self.run_summary(["BAMOK", "FQOK", "BAMGONE"],
                                   types={"BAMOK": "bam", "BAMGONE": "cram"})
        self.assertEqual(rows["BAMOK"]["terminal_state"], "completed_profiled")
        self.assertEqual(rows["FQOK"]["terminal_state"], "completed_profiled")
        self.assertEqual(rows["BAMGONE"]["terminal_state"], "failed_at_extraction")

    def test_stages_not_requested_are_not_expected(self):
        self.humann = []
        self.complete("S1", mpa=False, humann=False)
        rows, _ = self.run_summary(["S1"], stages="read_qc,host_depletion,classification")
        self.assertEqual(rows["S1"]["terminal_state"], "completed_profiled")
        self.assertEqual((rows["S1"]["metaphlan_status"], rows["S1"]["humann_status"]),
                         ("not_requested", "not_requested"))

    def test_host_depletion_only_run(self):
        self.humann = []
        fastp(self.d, "S1", 1000); depletion(self.d, "S1", 900)
        fastp(self.d, "S2", 1000)
        rows, _ = self.run_summary(["S1", "S2"], stages="extraction,read_qc,host_depletion")
        self.assertEqual(rows["S1"]["terminal_state"], "completed")
        self.assertEqual(rows["S2"]["terminal_state"], "failed_at_host_depletion")

    def test_humann_failure_when_every_sample_failed_it(self):
        """With --stages, a stage that EVERY sample failed is still a failure (cannot be inferred)."""
        self.humann = []
        self.complete("S1", humann=False)
        self.complete("S2", humann=False)
        rows, _ = self.run_summary(["S1", "S2"])
        self.assertEqual({r["terminal_state"] for r in rows.values()}, {"failed_at_humann"})

    def test_low_as_survival_only_for_rescue_on_samples(self):
        self.humann = []
        self.complete("ON")
        write(self.d / "ON.upstream_aligner.tsv",
              "sample\tupstream_aligner\tupstream_aligner_version\tupstream_bwa_min_score\tlow_as_mode\tlow_as_mode_reason\n"
              "ON\tbwa\t0.7\t0\ton\tauto_bwa_mem_T_below_30\n")
        write(self.d / "ON.extract_counts.tsv",
              "sample\tlibrary_primary_records\textracted_unmapped_records\tlow_as_reads_extracted\tlow_as_no_as_records\n"
              "ON\t1000\t100\t5\t0\n")
        rows, _ = self.run_summary(["ON"], types={"ON": "bam"})
        self.assertEqual(rows["ON"]["terminal_state"], "failed_at_low_as_survival")
        write(self.d / "ON.low_as_survival.tsv",
              "sample\tlow_as_reads_after_depletion\tlow_as_reads_classified\nON\t4\t3\n")
        rows, _ = self.run_summary(["ON"], types={"ON": "bam"})
        self.assertEqual(rows["ON"]["terminal_state"], "completed_profiled")

    def test_decontam_exclusion_and_failure_precedence(self):
        self.humann = []
        self.complete("KEEP")
        self.complete("LOWLIB")
        fastp(self.d, "FAILED", 1000)
        write(self.d / "x.sample_validation.tsv",
              "sample_id\tin_otu\tin_metadata\tsample_type\trole\tfinal_status\treason\n"
              "KEEP\tTRUE\tTRUE\tTumor\tpositive\tincluded\tincluded_positive\n"
              "LOWLIB\tTRUE\tTRUE\tTumor\tNA\texcluded\tinsufficient_library\n"
              "FAILED\tFALSE\tTRUE\tTumor\tNA\texcluded\tmetadata_only\n")
        rows, _ = self.run_summary(["KEEP", "LOWLIB", "FAILED"])
        self.assertEqual(rows["KEEP"]["terminal_state"], "completed_profiled")
        self.assertEqual(rows["LOWLIB"]["terminal_state"], "excluded_decontam_insufficient_library")
        # A failure upstream wins over decontam's "metadata_only": the failure is the cause.
        self.assertEqual(rows["FAILED"]["terminal_state"], "failed_at_host_depletion")

    def test_start_from_decontam_rows(self):
        self.humann = []
        write(self.d / "x.sample_validation.tsv",
              "sample_id\tin_otu\tin_metadata\tsample_type\trole\tfinal_status\treason\n"
              "A\tTRUE\tTRUE\tTumor\tpositive\tincluded\tincluded_positive\n"
              "B\tTRUE\tTRUE\tOther\tNA\texcluded\tunselected_sample_type\n")
        out = self.tmp / "o.tsv"
        subprocess.run([sys.executable, str(REPO / "scripts" / "build_qc_summary.py"), "--inputs",
                        str(self.tmp), "--output", str(out), "--stages", ""], check=True,
                       capture_output=True)
        with out.open() as handle:
            rows = {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(rows["A"]["terminal_state"], "completed")
        self.assertEqual(rows["B"]["terminal_state"], "excluded_decontam_unselected_sample_type")
        self.assertEqual(rows["A"]["metaphlan_status"], "NA")

    def test_exactly_one_state_and_known_vocabulary(self):
        self.humann = []
        self.complete("A")
        rows, _ = self.run_summary(["A", "B"])
        pattern = re.compile(r"^(completed(_profiled|_no_call)?|failed_at_(" + "|".join(QC.STAGE_ORDER)
                             + r")|excluded_[a-z0-9_]+)$")
        for row in rows.values():
            self.assertRegex(row["terminal_state"], pattern)

    def test_inferred_stages_without_the_flag(self):
        self.humann = []
        self.complete("A", humann=False)
        self.complete("B", mpa=False, humann=False)
        rows, _ = self.run_summary(["A", "B"], stages=None)
        self.assertEqual(rows["A"]["terminal_state"], "completed_profiled")
        self.assertEqual(rows["B"]["terminal_state"], "failed_at_metaphlan")

    def test_unknown_stage_is_an_error(self):
        with self.assertRaises(ValueError):
            QC.parse_stages("read_qc,assembly")


class RunReportOutcomes(unittest.TestCase):
    """RUN_REPORT gets counts per terminal state (no sample IDs)."""

    def test_sample_outcomes_counts_only(self):
        spec = importlib.util.spec_from_file_location("fp", REPO / "scripts" / "finalize_provenance.py")
        fp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fp)
        with tempfile.TemporaryDirectory() as tmp:
            qc = Path(tmp) / "cmp.qc.summary.tsv"
            qc.write_text("sample\tterminal_state\tmetaphlan_status\thumann_status\n"
                          "SECRET1\tcompleted_profiled\tran\tran\n"
                          "SECRET2\tfailed_at_classification\tNA\tNA\n"
                          "SECRET3\tfailed_at_classification\tNA\tNA\n")
            summary = fp.sample_outcomes({"qc_summary_dir": tmp})
            self.assertEqual(summary["by_state"], {"completed_profiled": 1, "failed_at_classification": 2})
            self.assertEqual(summary["failed"], 2)
            self.assertNotIn("SECRET", json.dumps(summary))
            text = "\n".join(fp.render_sample_outcomes(summary))
            self.assertIn("SAMPLE OUTCOMES", text)
            self.assertIn("failed_at_classification", text)
            self.assertIsNone(fp.sample_outcomes({"qc_summary_dir": str(Path(tmp) / "none")}))


if __name__ == "__main__":
    unittest.main()
