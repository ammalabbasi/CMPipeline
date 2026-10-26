#!/usr/bin/env python3
"""Database manifest and per-run database versions (audit X14; pksProfiler c08098d).

Synthetic database trees only (tiny files under a temp dir), plus a parse of the committed
conf/databases.tsv. The listing fingerprint is checked against the exact shell pipeline
pksProfiler's scripts/record_database_versions.sh uses, so the frozen values stay comparable.
"""

import csv
import hashlib
import importlib.util
import json
import os
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DV = load("database_versions", ROOT / "scripts" / "database_versions.py")
FP = load("finalize_provenance_dv", ROOT / "scripts" / "finalize_provenance.py")


def make_kraken(root: Path) -> Path:
    db = root / "kuniq"
    (db / "sub").mkdir(parents=True)
    (db / "database.kdb").write_bytes(b"K" * 5000)
    (db / "taxDB").write_text("1\t1\troot\tno rank\n")
    (db / "ktaxonomy.tsv").write_text("1\t|\t1\t|\troot\n")
    (db / "database.kdb.counts").write_text("1\t10\n")
    (db / "database50mers.kmer_distrib").write_text("mapped_taxid\tgenome_taxids:kmers_mapped:total_genome_kmers\n")
    (db / "database150mers.kmer_distrib").write_text("mapped_taxid\tother\n")
    (db / "sub" / "nested.txt").write_text("x")
    os.symlink(db / "taxDB", db / "link_to_taxDB")      # symlinks are not files of the database
    return db


def make_mmi(path: Path, seqs=(("chr1", 1000), ("chr2", 500)), tail=b"\x00" * 4096) -> bytes:
    """A minimal minimap2 index: magic, w k b n_seq flag, the dictionary, then 'buckets'."""
    head = b"MMI\x02" + struct.pack("<5I", 11, 21, 14, len(seqs), 0)
    for name, length in seqs:
        head += bytes([len(name)]) + name.encode() + struct.pack("<I", length)
    path.write_bytes(head + tail)
    return head


def manifest_for(root: Path, rows: list[dict]) -> Path:
    path = root / "databases.tsv"
    with path.open("w") as fh:
        fh.write("# comment lines are ignored\n")
        w = csv.DictWriter(fh, fieldnames=DV.MANIFEST_FIELDS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def freeze_row(name, param, path, select="*", key_method="", key_files="", **extra):
    r = DV.fingerprint(name, param, str(path), select, key_method, key_files)
    out = {f: "" for f in DV.MANIFEST_FIELDS}
    out.update({k: r[k] for k in ("size_bytes", "file_count", "fingerprint_method", "fingerprint",
                                  "newest_mtime", "key_sha256")})
    out.update(name=name, param=param, path=str(path), select=select, key_method=key_method,
               key_files=key_files, release=f"{name} release 1", frozen_on="2026-10-26", **extra)
    return out


class Fingerprints(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.db = make_kraken(self.tmp)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp)

    def test_listing_equals_pksprofiler_shell_method(self):
        """`find . -type f -printf '%P\\t%s\\n' | LC_ALL=C sort`, then printf '%s\\n' | sha256sum."""
        shell = ("listing=$(cd \"$1\" && find . -type f -printf '%P\\t%s\\n' | LC_ALL=C sort); "
                 "printf '%s\\n' \"$listing\" | sha256sum | cut -d' ' -f1")
        expected = subprocess.run(["bash", "-c", shell, "_", str(self.db)], capture_output=True,
                                  text=True, check=True).stdout.strip()
        row = DV.fingerprint("kraken_db", "kraken_db", str(self.db))
        self.assertEqual(row["fingerprint"], expected)
        self.assertEqual(row["file_count"], 7)              # the symlink is not counted
        self.assertEqual(row["fingerprint_method"], "sha256:listing")

    def test_single_large_file_is_listed_small_file_is_hashed(self):
        big = self.tmp / "big.mmi"
        make_mmi(big, tail=b"\x00" * 100)
        row = DV.fingerprint("x", "x", str(big))
        self.assertEqual(row["fingerprint_method"], "sha256:content")
        self.assertEqual(row["fingerprint"], hashlib.sha256(big.read_bytes()).hexdigest())
        old = DV.CONTENT_MAX
        DV.CONTENT_MAX = 10
        try:
            row = DV.fingerprint("x", "x", str(big))
        finally:
            DV.CONTENT_MAX = old
        self.assertEqual(row["fingerprint_method"], "sha256:listing")
        self.assertEqual(row["fingerprint"],
                         hashlib.sha256(f"big.mmi\t{big.stat().st_size}\n".encode()).hexdigest())

    def test_mmi_dictionary_key_reads_only_the_dictionary(self):
        mmi = self.tmp / "h.mmi"
        head = make_mmi(mmi)
        self.assertEqual(DV.mmi_dictionary_sha256(mmi), hashlib.sha256(head).hexdigest())
        make_mmi(mmi, tail=b"\x01" * 4096)                  # different buckets, same dictionary
        self.assertEqual(DV.mmi_dictionary_sha256(mmi), hashlib.sha256(head).hexdigest())
        make_mmi(mmi, seqs=(("chr1", 1000), ("chr2", 500), ("phiX", 5386)))
        self.assertNotEqual(DV.mmi_dictionary_sha256(mmi), hashlib.sha256(head).hexdigest())
        (self.tmp / "not.mmi").write_bytes(b"XXXX" + b"\x00" * 20)
        with self.assertRaises(ValueError):
            DV.mmi_dictionary_sha256(self.tmp / "not.mmi")

    def test_select_restricts_to_kmer_distrib(self):
        row = DV.fingerprint("bracken_kmer_distrib", "kraken_db", str(self.db), "database*mers.kmer_distrib",
                             "sha256:files", "database*mers.kmer_distrib")
        self.assertEqual(row["file_count"], 2)
        self.assertTrue(row["key_sha256"])

    def test_prefix_and_missing(self):
        (self.tmp / "idx.1.bt2l").write_text("a")
        (self.tmp / "idx.2.bt2l").write_text("bb")
        self.assertEqual(DV.fingerprint("x", "x", str(self.tmp / "idx"))["kind"], "prefix")
        self.assertEqual(DV.fingerprint("x", "x", str(self.tmp / "nothing"))["status"], "absent")
        self.assertEqual(DV.fingerprint("x", "x", None)["status"], "unset")


class Verdicts(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.db = make_kraken(self.tmp)
        self.mmi = self.tmp / "host.mmi"
        make_mmi(self.mmi)
        rows = [freeze_row("kraken_db", "kraken_db", self.db, key_method="sha256:files",
                           key_files="taxDB,ktaxonomy.tsv,database.kdb.counts", stage="classification"),
                freeze_row("bracken_kmer_distrib", "kraken_db", self.db, "database*mers.kmer_distrib",
                           "sha256:files", "database*mers.kmer_distrib", stage="classification"),
                freeze_row("hg38_db", "hg38_db", self.mmi, key_method="mmi:dictionary", stage="host_depletion")]
        self.manifest = manifest_for(self.tmp, rows)
        self.paths = {"kraken_db": str(self.db), "hg38_db": str(self.mmi), "metaphlan_db": None}

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp)

    def verdicts(self, stages=("host_depletion", "classification"), paths=None):
        out = DV.record(paths or self.paths, list(stages), self.manifest)
        return {r["name"]: r for r in out["rows"]}

    def test_frozen(self):
        rows = self.verdicts()
        for name in ("kraken_db", "bracken_kmer_distrib", "hg38_db"):
            self.assertEqual(rows[name]["verdict"], "frozen", rows[name])
        self.assertEqual(rows["kraken_db"]["release"], "kraken_db release 1")
        self.assertEqual(rows["metaphlan_db"]["verdict"], "unset")

    def test_added_file_changes_the_listing(self):
        (self.db / "database75mers.kmer_distrib").write_text("new")
        rows = self.verdicts()
        self.assertEqual(rows["kraken_db"]["verdict"], "CHANGED: listing")
        self.assertEqual(rows["bracken_kmer_distrib"]["verdict"], "CHANGED: listing")

    def test_same_size_edit_of_a_key_file_is_caught_by_the_key(self):
        text = (self.db / "taxDB").read_text()
        (self.db / "taxDB").write_text(text.replace("root", "ROOT"))      # same size
        rows = self.verdicts()
        self.assertEqual(rows["kraken_db"]["verdict"], "CHANGED: key")
        self.assertEqual(rows["bracken_kmer_distrib"]["verdict"], "frozen")

    def test_rebuilt_index_dictionary_is_caught(self):
        """A real .mmi is far above CMP_DB_HASH_MAX_BYTES, so only its listing and key are checked."""
        old = DV.CONTENT_MAX
        DV.CONTENT_MAX = 10
        try:
            manifest = manifest_for(self.tmp, [freeze_row("hg38_db", "hg38_db", self.mmi,
                                                          key_method="mmi:dictionary", stage="host_depletion")])
            make_mmi(self.mmi, seqs=(("chr1", 1000), ("chrX", 500)))      # same size, other contig
            rows = {r["name"]: r for r in DV.record({"hg38_db": str(self.mmi)}, ["host_depletion"], manifest)["rows"]}
        finally:
            DV.CONTENT_MAX = old
        self.assertEqual(rows["hg38_db"]["fingerprint_method"], "sha256:listing")
        self.assertEqual(rows["hg38_db"]["verdict"], "CHANGED: key")

    def test_touch_is_a_note_not_a_change(self):
        os.utime(self.db / "taxDB", (1, 4_000_000_000))
        row = self.verdicts()["kraken_db"]
        self.assertEqual(row["verdict"], "frozen")
        self.assertIn("mtime", row["note"])

    def test_copy_elsewhere_is_frozen_with_path_differs(self):
        import shutil
        copy = self.tmp / "copy"
        shutil.copytree(self.db, copy, symlinks=True)
        row = self.verdicts(paths={"kraken_db": str(copy)})["kraken_db"]
        self.assertEqual(row["verdict"], "frozen (path differs)")

    def test_stage_not_requested_is_unused_and_not_walked(self):
        rows = self.verdicts(stages=("host_depletion",))
        self.assertEqual(rows["kraken_db"]["verdict"], "unused")
        self.assertEqual(rows["kraken_db"]["fingerprint"], "")
        self.assertEqual(rows["hg38_db"]["verdict"], "frozen")

    def test_not_in_manifest_and_absent(self):
        other = self.tmp / "mpa"
        other.mkdir()
        (other / "mpa_latest").write_text("x")
        rows = self.verdicts(paths={"metaphlan_db": str(other), "kraken_db": str(self.tmp / "gone")})
        self.assertEqual(rows["metaphlan_db"]["verdict"], "not_in_manifest")
        self.assertEqual(rows["kraken_db"]["verdict"], "absent")

    def test_report_section(self):
        (self.db / "database75mers.kmer_distrib").write_text("new")
        lines = "\n".join(DV.render(DV.record(self.paths, ["host_depletion", "classification"], self.manifest)))
        self.assertIn("DATABASE VERSIONS", lines)
        self.assertIn("release  kraken_db release 1", lines)
        self.assertIn("NOT the frozen database", lines)


class FinaliseWritesIt(unittest.TestCase):
    def test_run_report_and_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db = make_kraken(tmp)
            tracedir = tmp / "pipeline_info"
            (tracedir / "provenance").mkdir(parents=True)
            record = {"status": "started", "run": {}, "environment": {},
                      "invocation": {"params": {"kraken_db": str(db), "hg38_db": "null"},
                                     "stages_requested": ["entry:beginning", "classification"]}}
            (tracedir / "provenance" / "s_r.json").write_text(json.dumps(record))
            argv = sys.argv
            sys.argv = ["finalize_provenance.py", "--tracedir", str(tracedir), "--session", "s_r",
                        "--success", "true"]
            try:
                self.assertEqual(FP.main(), 0)
            finally:
                sys.argv = argv
            written = json.loads((tracedir / "provenance" / "s_r.json").read_text())
            rows = {r["name"]: r for r in written["databases"]["rows"]}
            # Not the frozen TSCC database: a different tree under the same parameter.
            self.assertEqual(rows["kraken_db"]["verdict"], "CHANGED: listing")
            self.assertEqual(rows["validation_taxdb"]["verdict"], "unused")
            self.assertTrue(written["databases"]["manifest_sha256"])
            report = (tracedir / "RUN_REPORT.txt").read_text()
            self.assertIn("DATABASE VERSIONS", report)
            self.assertIn("kuniq_microbialdb 2023-08-08", report)   # the manifest's release label

    def test_never_fatal(self):
        out = FP.database_versions({"kraken_db": "/x"}, ["classification"], Path("/nonexistent/m.tsv"))
        self.assertEqual({r["name"]: r["verdict"] for r in out["rows"]}["kraken_db"], "absent")


class CommittedManifest(unittest.TestCase):
    def setUp(self):
        self.rows = DV.read_manifest(ROOT / "conf" / "databases.tsv")

    def test_every_tscc_database_is_frozen(self):
        tscc = (ROOT / "conf" / "tscc.config").read_text()
        for param in ("hg38_db", "t2t_phix_db", "kraken_db", "metaphlan_db", "humann3_nucleotide_db",
                      "humann3_protein_db"):
            with self.subTest(param=param):
                self.assertIn(f"params.{param}", tscc)
                self.assertIn(param, self.rows)
        self.assertIn("bracken_kmer_distrib", self.rows)

    def test_rows_are_complete_and_consistent(self):
        for name, row in self.rows.items():
            with self.subTest(name=name):
                for field in ("param", "stage", "release", "path", "fingerprint", "size_bytes", "frozen_on"):
                    self.assertTrue(row[field], field)
                self.assertIn(row["param"], DV.PARAM_STAGE)
                self.assertEqual(row["stage"], DV.PARAM_STAGE[row["param"]])
                self.assertRegex(row["fingerprint"], r"^[0-9a-f]{64}$")
                self.assertIn(row["fingerprint_method"], ("sha256:listing", "sha256:content"))
                if row["key_method"]:
                    self.assertIn(row["key_method"], ("sha256:files", "mmi:dictionary"))
                    self.assertRegex(row["key_sha256"], r"^[0-9a-f]{64}$")

    def test_shared_files_carry_pksprofilers_frozen_values(self):
        """Same files, same method: pksProfiler conf/databases.tsv (frozen 2026-10-14/22)."""
        self.assertEqual(self.rows["kraken_db"]["fingerprint"],
                         "fe43219f1104801e2a7e47a47209384def3fb312ef4e8a147c603c9135c9a8a6")
        self.assertEqual(self.rows["hg38_db"]["fingerprint"],
                         "92b1f812752bfe3097d302c0336cbb31135db729ab2361be5b30474fb43ab758")
        self.assertEqual(self.rows["t2t_phix_db"]["fingerprint"],
                         "8e2257962675aad5ef92adff8890084e51c77bc8e55e8b11b800d74ac3bfb196")

    def test_t2t_row_says_no_phix(self):
        self.assertIn("NO phiX", self.rows["t2t_phix_db"]["database"])

    def test_documented(self):
        outputs = (ROOT / "docs" / "outputs.md").read_text()
        self.assertIn("DATABASE VERSIONS", outputs)
        self.assertIn("conf/databases.tsv", outputs)
        self.assertIn("conf/databases.tsv", (ROOT / "docs" / "CHANGELOG.md").read_text())


if __name__ == "__main__":
    unittest.main()
