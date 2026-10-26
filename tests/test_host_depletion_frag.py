"""Host depletion is independent of input route (Ludmil E0 / audit X01, decision 37, 2026-10-06).

minimap2's -ax sr preset pairs ADJACENT mate-named reads (fragment mode). BAM/CRAM extraction emits
mates adjacent, paired FASTQ arrives as R1 block then R2 block, so the same reads were depleted
differently by route. Every host-depletion call must pass --frag=no unless --host_depletion_frag auto.
Behaviour on the real GRCh38/T2T indexes was measured with synthetic reads (SLURM 12425893).
"""
import re, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class HostDepletionFrag(unittest.TestCase):
    def test_default_is_no_and_validated(self):
        text = (ROOT / "main.nf").read_text()
        self.assertIn("params.host_depletion_frag = 'no'", text)
        self.assertIn("--host_depletion_frag must be no or auto", text)

    def test_every_depletion_call_carries_frag(self):
        mod = (ROOT / "Modules" / "map_reads.nf").read_text()
        calls = re.findall(r"minimap2 [^|\n]*-ax sr[^|\n]*", mod)
        self.assertEqual(len(calls), 3, calls)          # GRCh38, T2T-CHM13, pangenome loop
        for c in calls:
            self.assertIn("${frag}", c)
        self.assertIn("def frag = params.host_depletion_frag.toString() == 'no' ? '--frag=no' : ''", mod)

    def test_every_sr_alignment_in_the_pipeline_carries_frag(self):
        """2026-10-26: alignForValidation streamed interleaved mates into minimap2 -ax sr without
        --frag=no, the same route dependence as E0 (and pksProfiler's V8a bug, f134fb8)."""
        for nf in sorted((ROOT / "Modules").glob("*.nf")):
            code = "\n".join(l for l in nf.read_text().splitlines() if not l.lstrip().startswith(("*", "//", "/*", "#")))
            for c in re.findall(r"minimap2 [^|\n]*-ax sr[^|\n]*", code):
                self.assertIn("${frag}", c, f"{nf.name}: {c}")

    def test_documented(self):
        self.assertIn("`--host_depletion_frag`", (ROOT / "docs" / "parameters.md").read_text())


if __name__ == "__main__":
    unittest.main()
