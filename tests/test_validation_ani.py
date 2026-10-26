"""ANI curation of the validation reference and the group-level fallback (decision 36, 26 Oct 2026).

Synthetic only: ANI values are written by hand to mirror the Fusobacterium findings of 6 Oct 2026
(audits/2026-10-06_fuso_ani/report.md); the scorer test builds SAM records in code.
"""
import csv, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))
import validation_ani as VA   # noqa: E402
from test_alignment_validation import sam_line   # noqa: E402

G = 100_000


def pairs(spec):
    ani = {}
    for a, b, v in spec:
        ani[(a, b)] = ani[(b, a)] = v
    return ani


class Relabel(unittest.TestCase):
    """nuc1-3 true nucleatum; ani1-3 true animalis; mis = filed nucleatum, really animalis;
    odd1-2 = the unnamed cluster filed nucleatum (91.5% to nucleatum, < 93 to everything)."""

    def setUp(self):
        self.genomes = [("nuc1", "F_nucleatum"), ("nuc2", "F_nucleatum"), ("nuc3", "F_nucleatum"),
                        ("mis", "F_nucleatum"), ("odd1", "F_nucleatum"), ("odd2", "F_nucleatum"),
                        ("ani1", "F_animalis"), ("ani2", "F_animalis"), ("ani3", "F_animalis"),
                        ("sp1", "F_sp._X1")]
        spec = [("nuc1", "nuc2", 97.0), ("nuc1", "nuc3", 96.5), ("nuc2", "nuc3", 97.2),
                ("ani1", "ani2", 98.0), ("ani1", "ani3", 97.5), ("ani2", "ani3", 98.1),
                ("odd1", "odd2", 99.0)]
        for n in ("nuc1", "nuc2", "nuc3"):
            spec += [(n, "mis", 90.6), (n, "odd1", 91.5), (n, "odd2", 91.4)]
            spec += [(n, x, 91.0) for x in ("ani1", "ani2", "ani3")]
        for x in ("ani1", "ani2", "ani3"):
            spec += [(x, "mis", 97.7), (x, "odd1", 92.0), (x, "odd2", 92.1), (x, "sp1", 96.3)]
        self.ani = pairs(spec)

    def actions(self):
        return {r["accession"]: (r["action"], r["label"]) for r in VA.relabel(self.genomes, self.ani)}

    def test_mislabelled_genome_moves(self):
        self.assertEqual(self.actions()["mis"], ("relabel", "F_animalis"))

    def test_unnamed_cluster_stays_under_nucleatum(self):
        """Ammal, 26 Oct 2026: keep, do not split or drop."""
        got = self.actions()
        self.assertEqual(got["odd1"], ("divergent", "F_nucleatum"))
        self.assertEqual(got["odd2"], ("divergent", "F_nucleatum"))

    def test_seed_and_true_members_stay(self):
        got = self.actions()
        self.assertEqual(got["nuc1"], ("keep_seed", "F_nucleatum"))
        self.assertEqual(got["nuc2"], ("keep", "F_nucleatum"))
        self.assertEqual(got["ani2"], ("keep", "F_animalis"))

    def test_unnamed_label_is_a_source_not_a_target(self):
        got = self.actions()
        self.assertEqual(got["sp1"], ("relabel", "F_animalis"))   # an unnamed seed may move
        genomes = self.genomes + [("sp2", "F_sp._X1")]
        ani = dict(self.ani); ani.update(pairs([("sp2", x, 96.0) for x in ("ani1", "ani2", "ani3")]))
        got = {r["accession"]: (r["action"], r["label"]) for r in VA.relabel(genomes, ani)}
        self.assertEqual(got["sp2"], ("relabel", "F_animalis"))
        # nothing moves INTO an unnamed label
        self.assertTrue(all(lab != "F_sp._X1" for acc, (act, lab) in got.items() if act == "relabel"))

    def test_close_call_is_dropped(self):
        genomes = [("a1", "A_x"), ("a2", "A_x"), ("b1", "A_y"), ("b2", "A_y"), ("amb", "A_x")]
        ani = pairs([("a1", "a2", 97.0), ("b1", "b2", 97.0), ("amb", "a1", 95.5), ("amb", "a2", 95.5),
                     ("amb", "b1", 96.0), ("amb", "b2", 96.0)])
        got = {r["accession"]: r for r in VA.relabel(genomes, ani)}
        self.assertEqual((got["amb"]["action"], got["amb"]["label"]), ("drop_ambig", ""))


class Groups(unittest.TestCase):
    def test_92_joins_polymorphum_93_does_not(self):
        genomes = [("n", "F_nucleatum"), ("a", "F_animalis"), ("p", "F_polymorphum"), ("h", "F_hwasookii"),
                   ("o", "Other_x")]
        ani = pairs([("n", "a", 92.5), ("p", "h", 92.4), ("p", "n", 91.0), ("h", "n", 92.1)])
        taxid = {"F_nucleatum": "851", "F_animalis": "76859", "F_polymorphum": "76857", "F_hwasookii": "1583098"}
        g92 = VA.groups(genomes, ani, taxid, group_ani=92.0)
        self.assertEqual(set(g92), {"F_nucleatum", "F_animalis", "F_polymorphum", "F_hwasookii"})
        self.assertEqual(set(g92.values()), {"F_nucleatum_group"})     # named after the lowest taxid
        g93 = VA.groups(genomes, ani, taxid, group_ani=93.0)
        self.assertEqual(g93, {})


class GroupFallback(unittest.TestCase):
    """Two grouped species whose reads are mostly shared: species level cannot judge spread,
    the group can. A third grouped species with low identity must NOT be rescued."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        with open(f"{self.d}/contigs.tsv", "w") as f:
            f.write("contig\tlabel\taccession\tlength\tani_group\n")
            for lab, grp in (("F_ani", "F_nuc_group"), ("F_nuc", "F_nuc_group"), ("F_bad", "F_nuc_group"),
                             ("Solo_x", "")):
                f.write(f"{lab}|acc_{lab}|c1\t{lab}\tacc_{lab}\t{G}\t{grp}\n")
        with open(f"{self.d}/table.tsv", "w") as f:
            f.write("clade_name\tX\n")
            for sp in ("F ani", "F nuc", "F bad", "Solo x"):
                f.write(f"d__Bacteria|g__F|s__{sp}\t100\n")
        lines, n = [], 0

        def read(alns):
            nonlocal n
            n += 1
            for i, (ref, pos, AS, NM) in enumerate(alns):
                lines.append(sam_line(f"r{n}/1", 0 if i == 0 else 256, ref, pos, AS, NM))
        for i in range(10):          # 10 unique F_ani reads: too few to judge spread at species level
            read([("F_ani|acc_F_ani|c1", i * 9000, 200, 0)])
        for i in range(50):          # 50 reads shared F_ani / F_nuc, spread
            read([("F_ani|acc_F_ani|c1", 500 + i * 1900, 200, 0), ("F_nuc|acc_F_nuc|c1", 700 + i * 1900, 200, 0)])
        for i in range(60):          # F_bad: spread but 4% mismatches -> identity failure
            read([("F_bad|acc_F_bad|c1", i * 1600, 180, 4)])
        for i in range(10):          # ungrouped, too few
            read([("Solo_x|acc_Solo_x|c1", i * 9000, 200, 0)])
        with open(f"{self.d}/x.sam", "w") as f:
            f.write("@HD\tVN:1.6\n"); f.writelines(lines)

    def run_score(self, *extra):
        r = subprocess.run([sys.executable, str(ROOT / "scripts" / "validation_score.py"), "--sam", f"{self.d}/x.sam",
                            "--sample", "X", "--contigs", f"{self.d}/contigs.tsv",
                            "--species_table", f"{self.d}/table.tsv", "--out_dir", self.d, *extra],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        return {row["taxon"]: row for row in csv.DictReader(open(f"{self.d}/X.validation_stats.tsv"), delimiter="\t")}

    def test_shared_reads_validate_at_group_level(self):
        rows = self.run_score()
        self.assertEqual(rows["F_ani"]["status"], "validated_group")
        self.assertIn("unique_reads_for_bins", rows["F_ani"]["fail_reason"])   # why species level could not
        self.assertEqual(rows["F_ani"]["ani_group"], "F_nuc_group")
        self.assertEqual(rows["F_ani"]["group_status"], "validated")

    def test_identity_failure_is_never_rescued(self):
        rows = self.run_score()
        self.assertEqual(rows["F_bad"]["status"], "failed")
        self.assertIn("identity", rows["F_bad"]["fail_reason"])

    def test_ungrouped_taxon_unchanged(self):
        rows = self.run_score()
        self.assertEqual(rows["Solo_x"]["status"], "insufficient")
        self.assertEqual(rows["Solo_x"]["ani_group"], "NA")

    def test_shadow_rescued_only_on_request(self):
        rows = self.run_score()
        self.assertEqual(rows["F_nuc"]["status"], "failed")              # every read shared: a shadow
        self.assertIn("unique_fraction", rows["F_nuc"]["fail_reason"])
        rows = self.run_score("--group_rescue_shared")
        self.assertEqual(rows["F_nuc"]["status"], "validated_group")

    def test_reference_without_groups_behaves_as_before(self):
        with open(f"{self.d}/contigs.tsv") as f:
            body = [l.rsplit("\t", 1)[0] + "\n" for l in f]
        with open(f"{self.d}/contigs.tsv", "w") as f:
            f.writelines(body)
        rows = self.run_score()
        self.assertEqual(rows["F_ani"]["status"], "insufficient")


if __name__ == "__main__":
    unittest.main()
