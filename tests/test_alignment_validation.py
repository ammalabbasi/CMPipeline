"""Alignment validation scripts, on tiny synthetic alignments built in code (no real data).

Each scenario is one the 2026-09-23 simulation benchmark showed the scorer must get right.
"""
import csv, gzip, os, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / "scripts"
G = 100_000          # synthetic genome length


def sam_line(q, flag, ref, pos, AS=200, NM=0, seq="A" * 100):
    return "\t".join([q, str(flag), ref, str(pos + 1), "60", "100M", "*", "0", "0",
                      seq if not flag & 256 else "*", "*", f"AS:i:{AS}", f"NM:i:{NM}"]) + "\n"


class Scorer(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        with open(f"{self.d}/contigs.tsv", "w") as f:
            f.write("contig\tlabel\taccession\tlength\n")
            for lab in ("Gen_a", "Gen_b", "Gen_c", "Gen_d", "Gen_e", "Esc_coli", "Shi_flex"):
                f.write(f"{lab}|acc_{lab}|c1\t{lab}\tacc_{lab}\t{G}\n")
            f.write(f"Homo_sapiens|hs|chr1\tHomo_sapiens\ths\t1000000\n")
        with open(f"{self.d}/table.tsv", "w") as f:
            f.write("clade_name\tX\n")
            for sp in ("Gen a", "Gen b", "Gen c", "Gen d", "Gen e", "Esc coli", "Shi flex", "Not inref"):
                f.write(f"d__Bacteria|g__{sp.split()[0]}|s__{sp}\t100\n")
        with open(f"{self.d}/groups.tsv", "w") as f:
            f.write("group\tmember\nEscShi\tEsc coli\nEscShi\tShi flex\n")
        lines, n = [], 0
        def read(alns, seq="A" * 100):
            nonlocal n
            n += 1
            for i, (ref, pos, AS, NM) in enumerate(alns):
                lines.append(sam_line(f"r{n}/1", 0 if i == 0 else 256, ref, pos, AS, NM, seq))
        for i in range(60):   # present: spread across Gen_a
            read([("Gen_a|acc_Gen_a|c1", i * 1600, 200, 0)])
        for i in range(60):   # pile-up: all within 300 bp of Gen_b
            read([("Gen_b|acc_Gen_b|c1", 5000 + i * 5, 200, 0)])
        for i in range(60):   # relative: spread, 4% mismatches
            read([("Gen_c|acc_Gen_c|c1", i * 1600, 180, 4)])
        for i in range(60):   # shadow: every read equal-best on Gen_a and Gen_d
            read([("Gen_a|acc_Gen_a|c1", 800 + i * 1600, 200, 0), ("Gen_d|acc_Gen_d|c1", i * 1600, 200, 0)])
        for i in range(10):   # insufficient
            read([("Gen_e|acc_Gen_e|c1", i * 9000, 200, 0)])
        for i in range(60):   # E. coli / Shigella: shared, spread
            read([("Esc_coli|acc_Esc_coli|c1", i * 1600, 200, 0), ("Shi_flex|acc_Shi_flex|c1", i * 1600, 200, 0)])
        lines.append("\t".join(["unm", "4", "*", "0", "0", "*", "*", "0", "0", "A" * 100, "*"]) + "\n")
        with open(f"{self.d}/x.sam", "w") as f:
            f.write("@HD\tVN:1.6\n"); f.writelines(lines)

    def run_score(self, *extra):
        r = subprocess.run([sys.executable, str(S / "validation_score.py"), "--sam", f"{self.d}/x.sam",
                            "--sample", "X", "--contigs", f"{self.d}/contigs.tsv",
                            "--species_table", f"{self.d}/table.tsv", "--out_dir", self.d, *extra],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        return {row["taxon"]: row for row in csv.DictReader(open(f"{self.d}/X.validation_stats.tsv"), delimiter="\t")}

    def test_each_scenario_gets_the_right_status(self):
        rows = self.run_score()
        self.assertEqual(rows["Gen_a"]["status"], "validated")
        self.assertEqual(rows["Gen_b"]["status"], "failed"); self.assertIn("bins_ratio", rows["Gen_b"]["fail_reason"])
        self.assertEqual(rows["Gen_c"]["status"], "failed"); self.assertIn("identity", rows["Gen_c"]["fail_reason"])
        self.assertEqual(rows["Gen_d"]["status"], "failed"); self.assertIn("unique_fraction", rows["Gen_d"]["fail_reason"])
        self.assertEqual(rows["Gen_e"]["status"], "insufficient")
        self.assertEqual(rows["Not_inref"]["status"], "not_in_reference")
        # without grouping, E. coli and Shigella share every read: both fail as shadows
        self.assertEqual(rows["Esc_coli"]["status"], "failed")

    def test_group_scores_members_as_one_taxon(self):
        rows = self.run_score("--groups", f"{self.d}/groups.tsv")
        self.assertIn("EscShi", rows)
        self.assertEqual(rows["EscShi"]["status"], "validated")
        self.assertEqual(rows["EscShi"]["members"], "Esc_coli,Shi_flex")
        self.assertNotIn("Esc_coli", rows)

    def test_export_holds_only_validated_taxa(self):
        self.run_score("--export_reads")
        with gzip.open(f"{self.d}/X.validated_reads.fasta.gz", "rt") as fh:
            taxa = {l.strip().split("taxon=")[1] for l in fh if l.startswith(">")}
        self.assertEqual(taxa, {"Gen_a"})

    def test_batch_stream_matches_single_sample(self):
        single = self.run_score()
        body = [l for l in open(f"{self.d}/x.sam") if not l.startswith("@")]
        with open(f"{self.d}/multi.sam", "w") as f:
            for smp in ("P", "Q"):
                f.writelines(f"{smp}::{l}" for l in body)
        with open(f"{self.d}/table2.tsv", "w") as f:
            rows = open(f"{self.d}/table.tsv").read().splitlines()
            f.write("clade_name\tP\tQ\tZ\n")
            f.writelines(f"{r.split(chr(9))[0]}\t100\t100\t0\n" for r in rows[1:])
        r = subprocess.run([sys.executable, str(S / "validation_score.py"), "--sam", f"{self.d}/multi.sam",
                            "--samples", "P,Q,Z", "--contigs", f"{self.d}/contigs.tsv",
                            "--species_table", f"{self.d}/table2.tsv", "--out_dir", self.d],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        for smp in ("P", "Q"):
            got = {row["taxon"]: row["status"] for row in
                   csv.DictReader(open(f"{self.d}/{smp}.validation_stats.tsv"), delimiter="\t")}
            self.assertEqual(got, {t: v["status"] for t, v in single.items()})
        self.assertTrue(os.path.exists(f"{self.d}/Z.validation_stats.tsv"))   # no reads, still a file

    def test_non_contiguous_samples_are_refused(self):
        body = [l for l in open(f"{self.d}/x.sam") if not l.startswith("@")]
        with open(f"{self.d}/bad.sam", "w") as f:
            f.writelines(f"P::{l}" for l in body[:5]); f.writelines(f"Q::{l}" for l in body[:5])
            f.writelines(f"P::{l}" for l in body[5:10])
        with open(f"{self.d}/table2.tsv", "w") as f:
            f.write("clade_name\tP\tQ\nd__Bacteria|g__Gen|s__Gen a\t1\t1\n")
        r = subprocess.run([sys.executable, str(S / "validation_score.py"), "--sam", f"{self.d}/bad.sam",
                            "--samples", "P,Q", "--contigs", f"{self.d}/contigs.tsv",
                            "--species_table", f"{self.d}/table2.tsv", "--out_dir", self.d],
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not contiguous", r.stderr)


class AuditFixes(unittest.TestCase):
    """Audit 2026-09-23: B01 (underscore names), B02 (all-zero sample), B03 (spread shadows)."""

    def test_b01_bracken_underscore_names_match_taxdb(self):
        d = tempfile.mkdtemp()
        with open(f"{d}/taxDB", "w") as f:
            f.write("562\t561\tEscherichia coli\tspecies\n")
        with open(f"{d}/t.tsv", "w") as f:       # kreport2mpa writes s__Escherichia_coli
            f.write("clade_name\tX\nd__B|g__Escherichia|s__Escherichia_coli\t30\n")
        subprocess.run([sys.executable, str(S / "validation_species_list.py"), "--species_table", f"{d}/t.tsv",
                        "--taxdb", f"{d}/taxDB", "--out", f"{d}/o.tsv", "--key_out", f"{d}/k.txt"],
                       check=True, capture_output=True)
        row = next(csv.DictReader(open(f"{d}/o.tsv"), delimiter="\t"))
        self.assertEqual((row["label"], row["species"], row["taxid"]), ("Escherichia_coli", "Escherichia coli", "562"))

    def test_b01_ambiguous_label_gets_no_taxid(self):
        d = tempfile.mkdtemp()
        with open(f"{d}/taxDB", "w") as f:       # two species collide on one label
            f.write("111\t1\tFoo bar\tspecies\n222\t1\tFoo_bar\tspecies\n562\t561\tEscherichia coli\tspecies\n")
        with open(f"{d}/t.tsv", "w") as f:
            f.write("clade_name\tX\nd__B|g__Foo|s__Foo_bar\t30\nd__B|g__Escherichia|s__Escherichia_coli\t30\n")
        subprocess.run([sys.executable, str(S / "validation_species_list.py"), "--species_table", f"{d}/t.tsv",
                        "--taxdb", f"{d}/taxDB", "--out", f"{d}/o.tsv", "--key_out", f"{d}/k.txt"],
                       check=True, capture_output=True)
        rows = {r["label"]: r for r in csv.DictReader(open(f"{d}/o.tsv"), delimiter="\t")}
        self.assertEqual(rows["Foo_bar"]["taxid"], "NA")
        self.assertEqual(rows["Escherichia_coli"]["taxid"], "562")

    def test_b02_all_zero_sample_does_not_abort_the_mask(self):
        d = tempfile.mkdtemp()
        with open(f"{d}/table.tsv", "w") as f:
            f.write("clade_name\tX\tZ\nd__B|g__Gen|s__Gen_a\t10\t0\n")
        with open(f"{d}/X.validation_stats.tsv", "w") as f:
            f.write("sample\ttaxon\tstatus\nX\tGen_a\tvalidated\n")
        with open(f"{d}/Z.validation_stats.tsv", "w") as f:
            f.write("sample\ttaxon\tstatus\n")                    # header only: nothing to score
        r = subprocess.run([sys.executable, str(S / "apply_validation_mask.py"), "--species_table", f"{d}/table.tsv",
                            "--stats", f"{d}/X.validation_stats.tsv", f"{d}/Z.validation_stats.tsv", "--out_dir", d],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        row = next(csv.DictReader(open(f"{d}/validated.species.tsv"), delimiter="\t"))
        self.assertEqual((row["X"], row["Z"]), ("10", "0"))

    def test_b03_shadow_spread_over_many_genomes_fails_not_insufficient(self):
        d = tempfile.mkdtemp()
        with open(f"{d}/contigs.tsv", "w") as f:
            f.write("contig\tlabel\taccession\tlength\n")
            f.write(f"Real_sp|acc_r|c1\tReal_sp\tacc_r\t{G}\n")
            for g in range(60):
                f.write(f"Spa_sp|acc_{g}|c1\tSpa_sp\tacc_{g}\t{G}\n")
        with open(f"{d}/table.tsv", "w") as f:
            f.write("clade_name\tX\nd__B|g__Real|s__Real_sp\t100\nd__B|g__Spa|s__Spa_sp\t100\n")
        with open(f"{d}/x.sam", "w") as f:
            for i in range(200):     # every read equal-best on Real_sp and on one of 60 Spa_sp strains
                f.write(sam_line(f"r{i}", 0, "Real_sp|acc_r|c1", i * 450))
                f.write(sam_line(f"r{i}", 256, f"Spa_sp|acc_{i % 60}|c1", i * 450))
        r = subprocess.run([sys.executable, str(S / "validation_score.py"), "--sam", f"{d}/x.sam", "--sample", "X",
                            "--contigs", f"{d}/contigs.tsv", "--species_table", f"{d}/table.tsv", "--out_dir", d],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {x["taxon"]: x for x in csv.DictReader(open(f"{d}/X.validation_stats.tsv"), delimiter="\t")}
        self.assertEqual(rows["Spa_sp"]["status"], "failed")                  # was: insufficient
        self.assertIn("unique_fraction", rows["Spa_sp"]["fail_reason"])

    def test_b04_past_snapshot_is_refused_not_mislabelled(self):
        """Decision 26: a past RefSeq month cannot be re-downloaded, so the build must stop."""
        d = tempfile.mkdtemp()
        with open(f"{d}/species.tsv", "w") as f:
            f.write("label\tspecies\ttaxid\n")
        r = subprocess.run([sys.executable, str(S / "build_validation_ref.py"), "--species", f"{d}/species.tsv",
                            "--cache_dir", f"{d}/cache", "--out_dir", d, "--snapshot", "2000-01"],
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("past RefSeq snapshot cannot be re-downloaded", r.stderr)

    def test_b05_oversized_reference_is_refused_before_download(self):
        d = tempfile.mkdtemp(); cache = f"{d}/cache"; os.makedirs(cache)
        cols = ["assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid",
                "species_taxid", "organism_name", "infraspecific_name", "isolate", "version_status",
                "assembly_level", "release_type", "genome_rep", "seq_rel_date", "asm_name", "asm_submitter",
                "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq",
                "relation_to_type_material", "asm_not_live_date", "assembly_type", "group", "genome_size"]
        for g in ("bacteria", "archaea", "fungi", "protozoa"):   # = build_validation_ref.GROUPS: no network
            with open(f"{cache}/assembly_summary_{g}.txt", "w") as f:
                f.write("#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n")
                f.write("#" + "\t".join(cols) + "\n")
                if g == "bacteria":
                    row = dict.fromkeys(cols, "na")
                    row.update(assembly_accession="GCF_1.1", taxid="9", species_taxid="9", organism_name="Big one",
                               version_status="latest", assembly_level="Complete Genome",
                               refseq_category="reference genome", genome_size="60000000000",
                               ftp_path="https://example.invalid/GCF_1.1_x")
                    f.write("\t".join(row[c] for c in cols) + "\n")
        with open(f"{d}/species.tsv", "w") as f:
            f.write("label\tspecies\ttaxid\nBig_one\tBig one\t9\n")
        r = subprocess.run([sys.executable, str(S / "build_validation_ref.py"), "--species", f"{d}/species.tsv",
                            "--cache_dir", cache, "--out_dir", d, "--max_ref_gb", "40"],
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("above --validation_max_ref_gb", r.stderr)
        self.assertFalse(os.listdir(f"{cache}/genomes"))          # nothing was downloaded

    def test_unsupported_is_zeroed_decision_29(self):
        """Bracken >= 25 but < 10 supporting alignments: unsupported, and the mask zeroes it."""
        d = tempfile.mkdtemp()
        labs = {"Few_sp": (100, 3), "None_sp": (100, 0), "Some_sp": (100, 15), "Low_sp": (20, 0)}
        with open(f"{d}/contigs.tsv", "w") as f:
            f.write("contig\tlabel\taccession\tlength\n")
            for lab in labs:
                f.write(f"{lab}|acc_{lab}|c1\t{lab}\tacc_{lab}\t{G}\n")
        with open(f"{d}/table.tsv", "w") as f:
            f.write("clade_name\tX\n")
            for lab, (bracken, _) in labs.items():
                f.write(f"d__B|g__{lab.split('_')[0]}|s__{lab}\t{bracken}\n")
        with open(f"{d}/x.sam", "w") as f:
            for lab, (_, n) in labs.items():
                for i in range(n):
                    f.write(sam_line(f"{lab}_r{i}", 0, f"{lab}|acc_{lab}|c1", i * 5000))
        r = subprocess.run([sys.executable, str(S / "validation_score.py"), "--sam", f"{d}/x.sam", "--sample", "X",
                            "--contigs", f"{d}/contigs.tsv", "--species_table", f"{d}/table.tsv", "--out_dir", d],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {x["taxon"]: x["status"] for x in csv.DictReader(open(f"{d}/X.validation_stats.tsv"), delimiter="\t")}
        self.assertEqual(rows["Few_sp"], "unsupported")      # 3 of 100
        self.assertEqual(rows["None_sp"], "unsupported")     # 0 of 100: "none of it aligns"
        self.assertEqual(rows["Some_sp"], "insufficient")    # 15: some support, kept
        self.assertEqual(rows["Low_sp"], "insufficient")     # Bracken 20: too few to test, kept
        r = subprocess.run([sys.executable, str(S / "apply_validation_mask.py"), "--species_table", f"{d}/table.tsv",
                            "--stats", f"{d}/X.validation_stats.tsv", "--out_dir", d], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        masked = {x["clade_name"].rsplit("s__", 1)[1]: x["X"]
                  for x in csv.DictReader(open(f"{d}/validated.species.tsv"), delimiter="\t")}
        self.assertEqual((masked["Few_sp"], masked["None_sp"]), ("0", "0"))
        self.assertEqual((masked["Some_sp"], masked["Low_sp"]), ("100", "20"))

    def test_b09_human_row_is_never_scored_or_masked(self):
        """A Homo_sapiens row in a pass-through table: human is only the host competitor."""
        d = tempfile.mkdtemp()
        with open(f"{d}/contigs.tsv", "w") as f:
            f.write("contig\tlabel\taccession\tlength\n")
            f.write(f"Homo_sapiens|hs|chr1\tHomo_sapiens\ths\t{G}\n")
        with open(f"{d}/table.tsv", "w") as f:
            f.write("clade_name\tX\nd__Eukaryota|g__Homo|s__Homo_sapiens\t100\n")
        with open(f"{d}/x.sam", "w") as f:
            for i in range(10):      # few, piled-up reads: would 'fail' if human were scored
                f.write(sam_line(f"r{i}", 0, "Homo_sapiens|hs|chr1", 5000 + i))
        r = subprocess.run([sys.executable, str(S / "validation_score.py"), "--sam", f"{d}/x.sam", "--sample", "X",
                            "--contigs", f"{d}/contigs.tsv", "--species_table", f"{d}/table.tsv", "--out_dir", d],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {x["taxon"]: x for x in csv.DictReader(open(f"{d}/X.validation_stats.tsv"), delimiter="\t")}
        self.assertEqual(rows["Homo_sapiens"]["status"], "not_in_reference")


class Mask(unittest.TestCase):
    def test_only_failed_cells_are_zeroed_and_genus_rolls_up(self):
        d = tempfile.mkdtemp()
        with open(f"{d}/table.tsv", "w") as f:
            f.write("clade_name\tX\tY\n")
            f.write("d__B|g__Gen|s__Gen a\t10\t20\nd__B|g__Gen|s__Gen b\t5\t7\nd__B|g__Oth|s__Oth c\t3\t0\n")
        cols = "sample\ttaxon\tstatus\n"
        open(f"{d}/X.validation_stats.tsv", "w").write(cols + "X\tGen_a\tvalidated\nX\tGen_b\tfailed\nX\tOth_c\tinsufficient\n")
        open(f"{d}/Y.validation_stats.tsv", "w").write(cols + "Y\tGen_a\tfailed\nY\tGen_b\tvalidated\n")
        r = subprocess.run([sys.executable, str(S / "apply_validation_mask.py"), "--species_table", f"{d}/table.tsv",
                            "--stats", f"{d}/X.validation_stats.tsv", f"{d}/Y.validation_stats.tsv", "--out_dir", d],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        sp = {row["clade_name"]: row for row in csv.DictReader(open(f"{d}/validated.species.tsv"), delimiter="\t")}
        self.assertEqual((sp["d__B|g__Gen|s__Gen a"]["X"], sp["d__B|g__Gen|s__Gen a"]["Y"]), ("10", "0"))
        self.assertEqual((sp["d__B|g__Gen|s__Gen b"]["X"], sp["d__B|g__Gen|s__Gen b"]["Y"]), ("0", "7"))
        self.assertEqual(sp["d__B|g__Oth|s__Oth c"]["X"], "3")          # insufficient keeps its count
        ge = {row["clade_name"]: row for row in csv.DictReader(open(f"{d}/validated.genus.tsv"), delimiter="\t")}
        self.assertEqual((ge["d__B|g__Gen"]["X"], ge["d__B|g__Gen"]["Y"]), ("10", "7"))

    def test_missing_sample_stats_is_an_error(self):
        d = tempfile.mkdtemp()
        open(f"{d}/table.tsv", "w").write("clade_name\tX\tY\nd__B|g__Gen|s__Gen a\t1\t1\n")
        open(f"{d}/X.validation_stats.tsv", "w").write("sample\ttaxon\tstatus\nX\tGen_a\tvalidated\n")
        r = subprocess.run([sys.executable, str(S / "apply_validation_mask.py"), "--species_table", f"{d}/table.tsv",
                            "--stats", f"{d}/X.validation_stats.tsv", "--out_dir", d], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no validation stats", r.stderr)


class SpeciesList(unittest.TestCase):
    def test_support_cut_taxid_match_and_empty_list(self):
        d = tempfile.mkdtemp()
        open(f"{d}/taxDB", "w").write("562\t561\tEscherichia coli\tspecies\n851\t848\tFusobacterium nucleatum\tspecies\n")
        open(f"{d}/t.tsv", "w").write("clade_name\tX\tY\n"
                                      "d__B|g__Escherichia|s__Escherichia coli\t30\t0\n"
                                      "d__B|g__Fusobacterium|s__Fusobacterium nucleatum\t10\t20\n"
                                      "d__B|g__Mad|s__Mad eup\t99\t0\n")
        r = subprocess.run([sys.executable, str(S / "validation_species_list.py"), "--species_table", f"{d}/t.tsv",
                            "--taxdb", f"{d}/taxDB", "--out", f"{d}/o.tsv", "--key_out", f"{d}/k.txt"],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {r["species"]: r for r in csv.DictReader(open(f"{d}/o.tsv"), delimiter="\t")}
        self.assertEqual(set(rows), {"Escherichia coli", "Mad eup"})     # F. nucleatum max 20 < 25
        self.assertEqual(rows["Escherichia coli"]["taxid"], "562")
        self.assertEqual(rows["Mad eup"]["taxid"], "NA")
        self.assertEqual(len(open(f"{d}/k.txt").read().strip()), 16)
        r = subprocess.run([sys.executable, str(S / "validation_species_list.py"), "--species_table", f"{d}/t.tsv",
                            "--taxdb", f"{d}/taxDB", "--min_support", "1000", "--out", f"{d}/o.tsv",
                            "--key_out", f"{d}/k.txt"], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("nothing to validate", r.stderr)


class CompareLevels(unittest.TestCase):
    def test_agreement_labels(self):
        d = tempfile.mkdtemp()
        h = "taxon_id\tdetected_batch_count\tdetected_batches\tmin_batches\tfinal_contaminant\n"
        open(f"{d}/s.tsv", "w").write(h + "d__B|g__E|s__E coli\t2\t\t2\tTRUE\nd__B|g__P|s__P m\t0\t\t2\tFALSE\n")
        open(f"{d}/g.tsv", "w").write(h + "d__B|g__E\t2\t\t2\tTRUE\nd__B|g__B\t2\t\t2\tTRUE\n")
        subprocess.run([sys.executable, str(S / "compare_decontam_levels.py"), "--species_final", f"{d}/s.tsv",
                        "--genus_final", f"{d}/g.tsv", "--out", f"{d}/o.tsv"], check=True, capture_output=True)
        got = {r["taxon"]: r["agreement"] for r in csv.DictReader(open(f"{d}/o.tsv"), delimiter="\t")}
        self.assertEqual(got, {"d__B|g__E|s__E coli": "agree", "d__B|g__E": "agree", "d__B|g__B": "genus_only"})


class Wiring(unittest.TestCase):
    def test_module_and_params(self):
        text = (ROOT / "main.nf").read_text()
        self.assertIn("params.run_alignment_validation = false", text)
        self.assertIn("params.validation_min_identity = 0.98", text)
        self.assertIn("include { validationSpeciesList; buildValidationRef; alignForValidation;", text)
        self.assertIn("tuple('validated_species', t, meta, lib, 'species')", text)
        mod = (ROOT / "Modules" / "alignment_validation.nf").read_text()
        self.assertIn('grep -c "loaded/built the index"', mod)       # one-part index guard
        self.assertIn('storeDir { "${params.validation_cache_dir}/ref-${ref_key}" }', mod)


if __name__ == "__main__":
    unittest.main()
