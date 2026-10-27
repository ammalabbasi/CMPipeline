"""Focused behavioral tests for the June-Young integration corrections."""
import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from count_tables import read_counts
import test_decontam_zero_total as zero_tests


class CountAndMask(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, text):
        p = self.d / name
        p.write_text(text)
        return p

    def run_script(self, name, *args):
        return subprocess.run([sys.executable, str(ROOT / "scripts" / name), *map(str, args)],
                              cwd=self.d, capture_output=True, text=True)

    def test_library_totals_exclude_ancestor_rows(self):
        genus = self.write("g.tsv", "#Classification\tX.bracken.G.krakenreport.txt\nd__Bacteria\t200\nd__Bacteria|g__G\t200\n")
        species = self.write("s.tsv", "#Classification\tX.bracken.S.krakenreport.txt\nd__Bacteria\t200\nd__Bacteria|g__G\t200\nd__Bacteria|g__G|s__G_a\t200\n")
        r = self.run_script("prepare_bracken_tables.py", "--genus", genus, "--species", species,
                            "--out_genus", "g.out", "--out_species", "s.out", "--out_library", "lib.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(self.d / "lib.tsv") as f:
            row = next(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(row, {"sample": "X", "bracken_genus_total": "200", "bracken_species_total": "200"})
        self.assertEqual(read_counts(self.d / "s.out", "s")[0], ["X"])

    def test_invalid_count_schema_fails(self):
        for header, body in [
            ("X\tX", "1\t2"), ("X", ""), ("X", "NaN"), ("X", "Inf"),
            ("X", "-1"), ("X", "1.5"), ("X", "1\t2"), (" X", "1")
        ]:
            with self.subTest(header=header, body=body):
                p = self.write("bad.tsv", "clade_name\t" + header + "\nd__B|g__G|s__G_a\t" + body + "\n")
                with self.assertRaises(ValueError):
                    read_counts(p, "s")
        p = self.write("dup.tsv", "clade_name\tX\nd__B|g__G|s__G_a\t1\nd__B|g__G|s__G_a\t2\n")
        with self.assertRaises(ValueError):
            read_counts(p, "s")

    def test_canonical_filename_like_id_is_not_shortened(self):
        p = self.write("canonical.tsv", "clade_name\tX.bracken.S.krakenreport.txt\nd__B|g__G|s__G_a\t1\n")
        self.assertEqual(read_counts(p, "s")[0], ["X.bracken.S.krakenreport.txt"])

    def test_missing_or_invalid_positive_status_fails(self):
        table = self.write("table.tsv", "clade_name\tX\nd__B|g__G|s__G_a\t30\n")
        for body in ["", "X\tG_a\tbogus\n", "Y\tG_a\tvalidated\n",
                     "X\tG_a\tvalidated\nX\tG_a\tfailed\n"]:
            with self.subTest(body=body):
                stats = self.write("X.validation_stats.tsv", "sample\ttaxon\tstatus\n" + body)
                r = self.run_script("apply_validation_mask.py", "--species_table", table, "--stats", stats)
                self.assertNotEqual(r.returncode, 0, r.stdout)

    def test_legitimate_zero_and_insufficient_are_preserved(self):
        table = self.write("table.tsv", "clade_name\tX\tZ\nd__B|g__G|s__G_a\t12\t0\n")
        x = self.write("X.validation_stats.tsv", "sample\ttaxon\tstatus\nX\tG_a\tinsufficient\n")
        z = self.write("Z.validation_stats.tsv", "sample\ttaxon\tstatus\n")
        r = self.run_script("apply_validation_mask.py", "--species_table", table, "--stats", x, z)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(read_counts(self.d / "validated.species.tsv", "s")[1]["d__B|g__G|s__G_a"], {"X": 12, "Z": 0})

    def test_genus_sums_only_surviving_species(self):
        table = self.write("decontam.csv", "clade_name,X,Z\nd__B|g__G|s__G_a,10,0\nd__B|g__G|s__G_b,7,2\n")
        r = self.run_script("aggregate_decontaminated_species.py", "--input", table, "--output", "genus.csv",
                            "--excluded", "excluded.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(read_counts(self.d / "genus.csv", "g", delimiter=",")[1], {"d__B|g__G": {"X": 17, "Z": 2}})
        with open(self.d / "excluded.tsv") as f:
            self.assertEqual(list(csv.DictReader(f, delimiter="\t")), [])

    def test_species_without_genus_parent_is_listed_not_summed(self):
        # Lineages as they occur in real Bracken species tables.
        table = self.write("decontam.csv", "clade_name,X,Z\n"
                           "d__B|g__G|s__G_a,10,0\n"
                           "d__Bacteria|s__unidentified_bacterial_endosymbiont,4,0\n"
                           "d__B|p__P|c__C|o__O|f__F|s__Coxiella-like_endosymbiont,0,3\n"
                           "s__synthetic_construct,5,6\n"
                           "d__B|g__H|x__odd|s__H_a,2,2\n")
        r = self.run_script("aggregate_decontaminated_species.py", "--input", table, "--output", "genus.csv",
                            "--excluded", "excluded.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(read_counts(self.d / "genus.csv", "g", delimiter=",")[1], {"d__B|g__G": {"X": 10, "Z": 0}})
        with open(self.d / "excluded.tsv") as f:
            rows = {r["taxon"]: r for r in csv.DictReader(f, delimiter="\t")}
        self.assertEqual(set(rows), {"d__Bacteria|s__unidentified_bacterial_endosymbiont",
                                     "d__B|p__P|c__C|o__O|f__F|s__Coxiella-like_endosymbiont",
                                     "s__synthetic_construct", "d__B|g__H|x__odd|s__H_a"})
        self.assertEqual(rows["s__synthetic_construct"]["reason"], "no genus in lineage")
        self.assertEqual(rows["s__synthetic_construct"]["total_count"], "11")
        self.assertEqual(rows["s__synthetic_construct"]["samples_present"], "2")
        self.assertEqual(rows["d__B|g__H|x__odd|s__H_a"]["reason"], "genus is not the direct parent")

    @unittest.skipUnless(all(importlib.util.find_spec(m) for m in ("pandas", "seaborn", "adjustText")),
                         "consensus env (pandas/seaborn/adjustText) not available")
    def test_consensus_counts_species_without_genus(self):
        genus = self.write("g.tsv", "clade_name\tX\tZ\nd__B|g__G\t10\t2\n")
        species = self.write("s.tsv", "clade_name\tX\tZ\nd__B|g__G|s__G_a\t10\t2\n"
                             "d__Bacteria|s__unidentified_bacterial_endosymbiont\t4\t0\n"
                             "s__synthetic_construct\t5\t6\n")
        r = self.run_script("compute_consensus_taxa.py", "--bracken_genus", genus, "--bracken_species", species,
                            "--output", "plot.pdf", "--output_common_genus", "cg.tsv",
                            "--output_common_species", "cs.tsv", "--output_library_sizes", "lib.tsv",
                            "--output_status", "status.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(self.d / "status.tsv") as f:
            status = {row["key"]: row["value"] for row in csv.DictReader(f, delimiter="\t")}
        self.assertEqual(status["species_without_genus_dropped"], "2")
        self.assertEqual(status["species_rows_out"], "1")

    def test_genus_rollup_with_no_genus_at_all_fails(self):
        table = self.write("decontam.csv", "clade_name,X\ns__synthetic_construct,5\n")
        r = self.run_script("aggregate_decontaminated_species.py", "--input", table, "--output", "genus.csv",
                            "--excluded", "excluded.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("genus table would be empty", r.stderr)

    def test_root_count_missing_invalid_duplicate_and_true_zero(self):
        text = (ROOT / "Modules/Bracken.nf").read_text()
        shell = text[text.index("    MICROBIAL_READS="):text.index('    echo "Bracken microbial reads')]
        shell = shell.replace("\\$", "$")
        for body, ok in [
            ("#header\n", False),
            ("0\tbad\t0\t0\t0\t0\t1\tno rank\troot\n", False),
            ("0\t0\t0\t0\t0\t0\t1\tno rank\troot\n", True),
            ("0\t0\t0\t0\t0\t0\t1\tno rank\troot\n" * 2, False)
        ]:
            with self.subTest(body=body):
                report = self.write("kraken.tsv", body)
                r = subprocess.run(["bash", "-eu", "-c", shell], env=dict(os.environ, REPORT=str(report)),
                                   capture_output=True, text=True)
                self.assertEqual(r.returncode == 0, ok, r.stderr)

    def test_gate_is_microbial_reads(self):
        """Decision 1 (27 Oct 2026): root minus human (9606) and synthetic construct (32630)."""
        text = (ROOT / "Modules/Bracken.nf").read_text()
        shell = text[text.index("    MICROBIAL_READS="):text.index('    echo "Bracken microbial reads')]
        shell = shell.replace("\\$", "$") + '\necho "$MICROBIAL_READS"'
        body = ("1\t5000000\t0\t0\t0\t0\t1\tno rank\troot\n"
                "1\t1000\t0\t0\t0\t0\t9606\tspecies\t      Homo sapiens\n"
                "1\t4950000\t0\t0\t0\t0\t32630\tspecies\t    synthetic construct\n"
                "1\t49000\t0\t0\t0\t0\t2\tsuperkingdom\t  Bacteria\n")
        report = self.write("kraken.tsv", body)
        r = subprocess.run(["bash", "-eu", "-c", shell], env=dict(os.environ, REPORT=str(report)),
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout.strip().splitlines()[-1], "49000")

    def test_fastp_count_missing_invalid_and_true_zero(self):
        text = (ROOT / "Modules/filter_reads.nf").read_text()
        shell = text[text.index("    awk '"):text.rindex('    """')]
        shell = shell.replace("${sampleID}", "X").replace("\\$", "$")
        for value, ok in [("", False), ('"total_reads": "bad",', False),
                          ('"total_reads": -1,', False), ('"total_reads": 0,', True),
                          ('"total_reads": 123,', True)]:
            with self.subTest(value=value):
                self.write("X.fastp.json", '{\n"after_filtering": {\n' + value + '\n}\n}\n')
                r = subprocess.run(["bash", "-eu", "-c", shell], cwd=self.d,
                                   env=dict(os.environ, COUNT=str(self.d / "count.txt")), capture_output=True, text=True)
                self.assertEqual(r.returncode == 0, ok, r.stderr)


class StrictLibrary(unittest.TestCase):
    setUp = zero_tests.ZeroTotalBehaviour.setUp
    tearDown = zero_tests.ZeroTotalBehaviour.tearDown
    run_r = zero_tests.ZeroTotalBehaviour.run_r
    def test_bad_library_counts_and_ids_stop(self):
        for body in ["S1\tInf\n", "S1\t-1\n", "S1\t1.5\n", "S1\tbad\n",
                     "S1\t5000\nS1\t5000\n", " S1\t5000\n"]:
            with self.subTest(body=body):
                Path(self.dir, "bad_lib.tsv").write_text("sample\tbracken_genus_total\n" + body)
                r = self.run_r("--library_table", "bad_lib.tsv")
                self.assertNotEqual(r.returncode, 0)
                self.assertIn("ERROR", r.stderr)

    def test_all_uncallable_has_explicit_state(self):
        ids = [f"S{i}" for i in range(1, 13)]
        Path(self.dir, "otu.tsv").write_text("clade_name\t" + "\t".join(ids) +
            "\ng__Singleton\t5\t" + "\t".join(["0"] * 11) + "\n")
        r = self.run_r()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("insufficient_evidence", Path(self.dir, "t.run_manifest.tsv").read_text())
        self.assertIn("retained_not_callable", Path(self.dir, "t.filter_report.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
