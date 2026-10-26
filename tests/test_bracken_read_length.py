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


class GroovyStringEscapes(unittest.TestCase):
    """Every backslash in a module's triple-quoted script must be a Groovy string escape.

    Bracken.nf's kmer_distrib sed was written with single backslashes before ( ) . and 1. Nextflow's
    v1 parser refused to compile main.nf; the v2 parser (26.04) compiled it and shifted every later
    interpolation of the Bracken script onto the wrong value, so every Bracken task failed. Neither
    -preview nor the emulation above (which un-escapes as Groovy would) could see it.
    """
    ALLOWED = set('\\$ntr"\'\nbfu0/')

    def offending(self, text):
        found = []
        for match in re.finditer(r'"""(.*?)"""', text, re.S):
            body, i = match.group(1), 0
            while i < len(body):
                if body[i] == "\\":
                    nxt = body[i + 1] if i + 1 < len(body) else ""
                    if nxt not in self.ALLOWED:
                        found.append((text[:match.start(1) + i].count("\n") + 1, body[max(0, i - 25):i + 15]))
                    i += 2
                else:
                    i += 1
        return found

    def test_module_scripts_use_valid_escapes(self):
        for path in sorted((ROOT / "Modules").glob("*.nf")) + [ROOT / "main.nf"]:
            with self.subTest(path=path.name):
                self.assertEqual(self.offending(path.read_text()), [])

    def test_the_scan_catches_the_old_line(self):
        old = 'x = """\n    sed -n \'s/^database\\([0-9]*\\)mers\\.k\\$/\\1/p\'\n"""'
        self.assertEqual(len(self.offending(old)), 4)


if __name__ == "__main__":
    unittest.main()
