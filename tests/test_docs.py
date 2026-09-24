"""The docs (README.md + docs/**/*.md) document every parameter main.nf defines, and have not
drifted back (audit A17/C18).

Since 2026-09-24 the README is a short front page (pksProfiler layout) and the detail lives in
docs/running/*.md, docs/hpc.md, docs/outputs.md and docs/parameters.md; the checks read all of them.
Parses every `params.<name> =` assignment in main.nf and asserts the docs mention `--<name>`.
A parameter added to main.nf without a README line fails here, which is how today's decontam and
batch-correction parameters went undocumented. Also pins a few statements the audit found wrong.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN_NF = ROOT / "main.nf"
README = ROOT / "README.md"
DOCS = [README, *sorted((ROOT / "docs").rglob("*.md"))]


def docs_text():
    return "\n".join(p.read_text() for p in DOCS)

# Parameters deliberately NOT documented in the README. Each entry needs a reason.
# Empty on purpose (audit A17): every param in main.nf, including the internal script/env/dir
# paths, is listed under the README's "All parameters" table, because each can be overridden on
# the command line and a user who meets one in a run record should be able to look it up.
UNDOCUMENTED_OK = {}

ASSIGNMENT = re.compile(r"^\s*params\.([A-Za-z_][A-Za-z0-9_]*)\s*=(?!=)", re.MULTILINE)


def main_nf_params():
    return sorted(set(ASSIGNMENT.findall(MAIN_NF.read_text())))


class ReadmeDocumentsEveryParam(unittest.TestCase):
    def test_parser_finds_params(self):
        # Guard the regex itself: an empty or tiny list would make the check below vacuous.
        found = main_nf_params()
        self.assertGreater(len(found), 50)
        for name in ("decontam_zero_total_policy", "batch_corr_method", "validation_min_support",
                     "hg38_db", "outdir"):
            self.assertIn(name, found)

    def test_every_param_is_in_readme(self):
        readme = docs_text()
        missing = [n for n in main_nf_params()
                   if n not in UNDOCUMENTED_OK
                   and not re.search(r"--" + re.escape(n) + r"(?![A-Za-z0-9_])", readme)]
        self.assertEqual(missing, [], "params in main.nf with no --<name> in README.md or docs/")

    def test_parameter_reference_is_complete(self):
        """docs/parameters.md is the one table listing every parameter."""
        ref = (ROOT / "docs" / "parameters.md").read_text()
        missing = [n for n in main_nf_params() if f"`--{n}`" not in ref]
        self.assertEqual(missing, [], "params missing from docs/parameters.md")

    def test_readme_links_resolve(self):
        import re as _re
        for p in DOCS:
            for target in _re.findall(r"\]\(([^)#:]+)\)", p.read_text()):
                with self.subTest(doc=p.name, link=target):
                    self.assertTrue((p.parent / target).exists(), f"{p.name} links to missing {target}")
        self.assertIn('src="workflow_logo/v0.2.png"', README.read_text())
        self.assertTrue((ROOT / "workflow_logo" / "v0.2.png").exists())

    def test_allow_list_is_not_stale(self):
        stale = sorted(set(UNDOCUMENTED_OK) - set(main_nf_params()))
        self.assertEqual(stale, [], "UNDOCUMENTED_OK names params main.nf no longer defines")


class ReadmeMatchesCode(unittest.TestCase):
    """Statements audit A17/C18 found contradicting the code."""

    @classmethod
    def setUpClass(cls):
        cls.readme = docs_text()
        cls.main = MAIN_NF.read_text()

    def test_unevaluable_default_is_skip(self):
        self.assertRegex(self.main, r'params\.decontam_unevaluable_policy\s*=\s*"skip"')
        self.assertRegex(self.readme, r"`--decontam_unevaluable_policy`\s*\|\s*\**`skip`")
        self.assertNotIn("unevaluable now stops the run", self.readme)

    def test_humann_skips_below_gate_samples(self):
        self.assertNotIn("computes its own", self.readme)
        self.assertIn("HUMAnN3 runs only on samples that cleared `--consensus_min_reads`", self.readme)

    def test_batch_corr_method_is_tune_or_vanilla(self):
        self.assertIn("if (!(params.batch_corr_method in ['tune', 'vanilla']))", self.main)
        self.assertNotRegex(self.readme.lower(), r"combat_seq|--batch_corr_method combat")

    def test_no_stale_version_or_gate_or_fork(self):
        for stale in ("25.10.4", "minus unclassified", "wdl2459", "nf-batch-correction-env"):
            self.assertNotIn(stale, self.readme)

    def test_qc_column_count_matches_script(self):
        import sys
        sys.path.insert(0, str(ROOT / "scripts"))
        import build_qc_summary
        n = len(build_qc_summary.COLUMNS)
        counts = {int(x) for x in re.findall(r"(\d+)[ -]columns?", self.readme)}
        self.assertEqual(counts, {n})
        for col in ("library_primary_records", "extracted_unmapped_records", "pct_retained_of_library"):
            self.assertIn(col, build_qc_summary.COLUMNS)
            self.assertIn(col, self.readme)

    def test_consensus_status_documented(self):
        self.assertIn("consensus_status.tsv", self.readme)
        self.assertIn("bracken_passthrough", self.readme)


class ConqurPinIsConsistent(unittest.TestCase):
    """The ConQuR fork and commit named in the yml comment and installer match params.conqur_sha."""

    def test_same_fork_and_sha(self):
        sha = re.search(r'params\.conqur_sha\s*=\s*"([0-9a-f]{40})"', MAIN_NF.read_text()).group(1)
        yml = (ROOT / "conda_envs" / "batch_correction_env.yml").read_text()
        installer = (ROOT / "scripts" / "install_conqur.R").read_text()
        self.assertIn(f'unset = "{sha}"', installer)
        self.assertIn('unset = "ivartb/ConQuR_par"', installer)
        self.assertIn("ivartb/ConQuR_par@" + sha[:8], yml)
        self.assertNotIn("wdl2459", yml)


if __name__ == "__main__":
    unittest.main()
