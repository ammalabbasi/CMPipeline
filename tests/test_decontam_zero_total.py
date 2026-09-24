"""Zero-total samples in decontamination (added 2026-09-23).

decontam::isContaminant() drops every sample with zero total counts. On the balanced CRC run
that removed all three batch-1 normals (none of the 6 consensus genera was present in them)
and with them the batch's controls. These tests pin the fix: all-zero samples are kept as
evidence of absence unless their library was empty upstream.

Static checks only (no R needed); the behavioural test is test_zero_total_behaviour, which
runs the R script when a decontam environment is available and is skipped otherwise.
"""
import os, shutil, subprocess, tempfile, unittest, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class ZeroTotalWiring(unittest.TestCase):
    def test_policy_parameter_defaults_to_keep_and_is_validated(self):
        text = (ROOT / "main.nf").read_text()
        self.assertIn('params.decontam_zero_total_policy = "keep"', text)
        self.assertIn('params.decontam_unevaluable_policy = "skip"', text)
        self.assertIn("--decontam_unevaluable_policy must be one of: skip, error", text)
        self.assertIn("--decontam_zero_total_policy must be one of: keep, drop", text)

    def test_consensus_emits_library_sizes_and_decontam_receives_them(self):
        module = (ROOT / "Modules" / "preprocess_taxa.nf").read_text()
        self.assertIn('path "bracken_library_sizes.tsv"', module)
        self.assertIn("--output_library_sizes bracken_library_sizes.tsv", module)
        self.assertIn("output_library_sizes", (ROOT / "scripts" / "compute_consensus_taxa.py").read_text())
        text = (ROOT / "main.nf").read_text()
        self.assertIn("consensus_taxa.out[4].set { CONSENSUS_LIBRARY_SIZES }", text)
        self.assertIn(".combine(CONSENSUS_LIBRARY_SIZES)", text)
        decon = (ROOT / "Modules" / "decontamination.nf").read_text()
        self.assertIn("path(library_sizes, stageAs:", decon)
        self.assertIn('--zero_total_policy "${params.decontam_zero_total_policy}"', decon)

    def test_r_script_uses_per_taxon_test_and_self_checks_it(self):
        r = (ROOT / "scripts" / "decontamination.R").read_text()
        self.assertIn('get0("isContaminantPrevalence", envir = asNamespace("decontam")', r)
        self.assertIn("self-check failed", r)
        self.assertIn("insufficient_library", r)
        self.assertIn("zero_total_input", r)


def _decontam_rscript():
    for env in glob.glob(str(ROOT / ".conda_cache" / "env-*")):
        if os.path.isdir(os.path.join(env, "lib", "R", "library", "decontam")):
            return os.path.join(env, "bin", "Rscript")
    return None


@unittest.skipUnless(_decontam_rscript(), "no decontam conda env built yet")
class ZeroTotalBehaviour(unittest.TestCase):
    """Synthetic 2-batch cohort, 3 tumours + 3 normals each, batch-1 normals all zero."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        S = [f"S{i}" for i in range(1, 13)]
        batch = {s: "batch1" if i < 6 else "batch2" for i, s in enumerate(S)}
        typ = {s: "Tumor" if i % 6 < 3 else "Normal" for i, s in enumerate(S)}
        vals = [[(7 * i + 11 * j) % 50 + 5 for i in range(12)] for j in range(6)]
        with open(f"{self.dir}/otu.tsv", "w") as f:
            f.write("clade_name\t" + "\t".join(S) + "\n")
            for j in range(6):
                row = ["0" if batch[s] == "batch1" and typ[s] == "Normal" else str(vals[j][i])
                       for i, s in enumerate(S)]
                f.write(f"g__G{j}\t" + "\t".join(row) + "\n")
        with open(f"{self.dir}/meta.tsv", "w") as f:
            f.write("patient\tsample_type\tshipment_batch\n")
            for s in S:
                f.write(f"{s}\t{typ[s]}\t{batch[s]}\n")
        def lib(name, sizes):
            with open(f"{self.dir}/{name}", "w") as f:
                f.write("sample\tbracken_genus_total\tbracken_species_total\n")
                for s in S:
                    f.write(f"{s}\t{sizes.get(s, 5000)}\t0\n")
        lib("lib_all.tsv", {})                 # every library comfortably above --min_library
        lib("lib.tsv", {"S4": 0})              # one truly empty library
        lib("lib_shallow.tsv", {"S5": 3})      # one near-empty library (3 reads)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def run_r(self, *extra):
        cmd = [_decontam_rscript(), str(ROOT / "scripts" / "decontamination.R"),
               "--otu_table", "otu.tsv", "--metadata", "meta.tsv", "--prefix", "t",
               "--start_mode", "decontam", "--environment_spec", str(ROOT / "conda_envs" / "decontam_env.yml"),
               "--nextflow_version", "test", "--taxon_column", "clade_name", "--taxonomic_rank", "genus",
               "--sample_id_column", "patient", "--batch_column", "shipment_batch",
               "--type_column", "sample_type", "--control_role", "surrogate_biological",
               "--positive_values", "Tumor", "--control_values", "Normal",
               "--unselected_type_policy", "error", "--threshold", "0.1", "--min_prevalence", "0.05",
               "--min_abundance", "5", "--min_batches", "2", "--min_total_per_batch", "5",
               "--min_positive_per_batch", "1", "--min_control_per_batch", "1", *extra]
        if "--library_table" not in extra and "--allow_unknown_library" not in extra and "--no_lib" not in extra:
            cmd += ["--library_table", "lib_all.tsv"]
        cmd = [c for c in cmd if c != "--no_lib"]
        return subprocess.run(cmd, cwd=self.dir, capture_output=True, text=True)

    def test_keep_retains_zero_controls(self):
        r = self.run_r()
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        summary = open(f"{self.dir}/t.decontamination_summary.txt").read()
        self.assertIn("Selected controls: 6", summary)
        self.assertIn("Eligible batches: 2", summary)

    def test_near_empty_library_is_not_evidence_of_absence(self):
        """A06: a 3-read library can't show presence; below --min_library it's excluded."""
        r = self.run_r("--library_table", "lib_shallow.tsv")
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        rows = [l.rstrip("\n").split("\t") for l in open(f"{self.dir}/t.sample_validation.tsv")]
        s5 = next(x for x in rows[1:] if x[0] == "S5")
        self.assertEqual(s5[rows[0].index("reason")], "insufficient_library")

    def test_all_zero_after_abundance_filter_gets_the_library_check(self):
        """Decision 24 (A06): a tiny library whose reads all fall below --min_abundance is all-zero
        in what the prevalence test sees, so it must face --min_library too."""
        rows = [l.rstrip("\n").split("\t") for l in open(f"{self.dir}/otu.tsv")]
        col = rows[0].index("S10")                           # a batch-2 control
        for r in rows[1:]:
            r[col] = "2"                                     # every count below --min_abundance 5
        with open(f"{self.dir}/otu.tsv", "w") as f:
            f.write("\n".join("\t".join(r) for r in rows) + "\n")
        with open(f"{self.dir}/lib_tiny.tsv", "w") as f:
            f.write("sample\tbracken_genus_total\tbracken_species_total\n")
            for s in [f"S{i}" for i in range(1, 13)]:
                f.write(f"{s}\t{12 if s == 'S10' else 5000}\t0\n")
        r = self.run_r("--library_table", "lib_tiny.tsv")
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        out = [l.rstrip("\n").split("\t") for l in open(f"{self.dir}/t.sample_validation.tsv")]
        s10 = next(x for x in out[1:] if x[0] == "S10")
        self.assertEqual(s10[out[0].index("reason")], "insufficient_library")

    def test_unknown_library_stops_unless_allowed(self):
        """A07: a zero-total sample with no library size must not be silently kept."""
        r = self.run_r("--no_lib")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("have no library size", r.stderr)
        r = self.run_r("--allow_unknown_library")
        self.assertEqual(r.returncode, 0, r.stderr[-500:])

    def test_taxon_table_as_library_source(self):
        """A10/A07: with consensus skipped, the unmasked table itself gives the library sizes;
        an all-zero sample there has library 0, so it is excluded as insufficient_library."""
        r = self.run_r("--library_table", "otu.tsv")
        rows = [l.rstrip("\n").split("\t") for l in open(f"{self.dir}/t.sample_validation.tsv")] \
            if os.path.exists(f"{self.dir}/t.sample_validation.tsv") else None
        # batch 1 loses its three normals (library 0), so decontam stops on batch eligibility,
        # naming the exclusions -- which is the correct outcome for truly empty libraries
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("3 zero-total sample(s) excluded", r.stderr)

    def test_contaminants_final_says_whether_a_taxon_was_callable(self):
        """A19: n_evaluable_batches and callable columns."""
        r = self.run_r()
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        head = open(f"{self.dir}/t.contaminants_final.tsv").readline().rstrip("\n").split("\t")
        self.assertIn("n_evaluable_batches", head)
        self.assertIn("callable", head)

    def test_drop_reproduces_the_old_failure_with_counts(self):
        r = self.run_r("--zero_total_policy", "drop")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("Eligible batches 1 below required 2", r.stderr)
        self.assertIn("3 zero-total sample(s) excluded", r.stderr)

    def add_singleton_taxon(self):
        """A taxon present in exactly one sample of batch 2 (S7): untestable there."""
        with open(f"{self.dir}/otu.tsv", "a") as f:
            f.write("g__Single\t" + "\t".join("40" if i == 6 else "0" for i in range(12)) + "\n")

    def test_unevaluable_taxon_is_recorded_not_fatal(self):
        self.add_singleton_taxon()
        r = self.run_r()
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        self.assertIn("unevaluable", r.stdout)
        summary = open(f"{self.dir}/t.decontamination_summary.txt").read()
        self.assertRegex(summary, r"Unevaluable taxon/batch pairs \(policy skip\): [1-9]")

    def test_unevaluable_policy_error_keeps_the_strict_behaviour(self):
        self.add_singleton_taxon()
        r = self.run_r("--unevaluable_policy", "error")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("Unevaluable decontam results", r.stderr)

    def test_empty_library_is_the_only_exclusion(self):
        r = self.run_r("--library_table", "lib.tsv")
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        rows = [l.rstrip("\n").split("\t") for l in open(f"{self.dir}/t.sample_validation.tsv")]
        head = rows[0]
        s4 = next(r for r in rows[1:] if r[0] == "S4")
        s5 = next(r for r in rows[1:] if r[0] == "S5")
        self.assertEqual(s4[head.index("reason")], "insufficient_library")
        self.assertEqual(s5[head.index("final_status")], "included")


if __name__ == "__main__":
    unittest.main()
