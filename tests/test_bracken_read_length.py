"""Bracken read length is chosen per sample (2026-10-26; pksProfiler c6599c8).

A fixed --bracken_read_length 50 was applied to ~150 bp TCGA WGS reads. `auto` (default) takes the
median length of the sample's reads and the database kmer_distrib length closest to it (ties ->
shorter); an integer still forces one length. The shell block is executed on synthetic reads.
"""
import os, re, subprocess, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MOD = (ROOT / "Modules" / "Bracken.nf").read_text()


def block():
    start = MOD.index("    MEDIAN_LEN=")
    end = MOD.index("    for lvl in G S; do")
    return MOD[start:end]


def run(lengths, db_lengths, setting="auto"):
    d = tempfile.mkdtemp()
    db = Path(d, "db"); db.mkdir()
    for L in db_lengths:
        (db / f"database{L}mers.kmer_distrib").write_text("")
    with open(Path(d, "S1.trimmed.fastq"), "w") as fh:
        for i, L in enumerate(lengths):
            fh.write(f"@r{i}\n{'A' * L}\n+\n{'I' * L}\n")
    sh = (block().replace("${params.kraken_db}", str(db)).replace("${SAMPLE_NAME}", "S1")
          .replace("${params.bracken_read_length}", str(setting)).replace("\\$", "$").replace("\\\\", "\\"))
    r = subprocess.run(["bash", "-euc", sh], cwd=d, capture_output=True, text=True)
    out = Path(d, "S1.bracken_read_length.tsv")
    rows = [l.split("\t") for l in out.read_text().splitlines()] if out.exists() else None
    return r, (dict(zip(rows[0], rows[1])) if rows else None)


class BrackenReadLength(unittest.TestCase):
    DB = [50, 75, 100, 125, 150, 200, 250, 300]

    def test_default_is_auto_and_validated(self):
        text = (ROOT / "main.nf").read_text()
        self.assertIn("params.bracken_read_length = 'auto'", text)
        self.assertIn("--bracken_read_length must be auto or a positive integer", text)
        self.assertIn("-r \\${BRACKEN_LEN}", MOD)
        self.assertNotIn("-r ${params.bracken_read_length}", MOD)

    def test_150bp_reads_use_150(self):
        r, row = run([151] * 50 + [120] * 10, self.DB)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(row["median_read_length"], "151")
        self.assertEqual(row["bracken_read_length"], "150")
        self.assertEqual(row["mode"], "auto")

    def test_tie_goes_to_shorter(self):
        _, row = run([175] * 9, self.DB)
        self.assertEqual(row["bracken_read_length"], "150")

    def test_short_reads_use_smallest_available(self):
        _, row = run([36] * 9, self.DB)
        self.assertEqual(row["bracken_read_length"], "50")

    def test_fixed_value_is_respected(self):
        _, row = run([151] * 9, self.DB, setting=100)
        self.assertEqual(row["bracken_read_length"], "100")
        self.assertEqual(row["mode"], "fixed")

    def test_no_distrib_fails(self):
        r, _ = run([151] * 9, [])
        self.assertNotEqual(r.returncode, 0)


if __name__ == "__main__":
    unittest.main()
