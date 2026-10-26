#!/usr/bin/env python3
"""The --t2t_phix_db index contains no phiX; nothing may say it removes phiX (port review #19).

human-GCA-phix-db.mmi is T2T-CHM13v2.0 alone (25 sequences, read from the index by pksProfiler on
2026-10-22; docs/host_depletion_references.md). The parameter and column names stay; the labels a
reader sees must say T2T-CHM13. Code and docs only -- the index is not read here.
"""
import importlib.util
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_cohort_report_labels", ROOT / "scripts" / "build_cohort_report.py")
CR = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CR)

# Places allowed to quote the old wording: the ones that explain it was wrong.
EXPLAINERS = {"docs/CHANGELOG.md", "docs/host_depletion_references.md", "tests/test_host_depletion_labels.py"}
CLAIM = re.compile(r"T2T\s*(\+|and)\s*phi\s*X|CHM13\s*\+\s*phi\s*X|Removed:[^\"\n]*phi\s*X", re.I)


class Labels(unittest.TestCase):
    def test_cohort_report_stage_label(self):
        labels = {key: label for key, label, _ in CR.STAGES}
        self.assertEqual(labels["pct_removed_t2t_phix"], "Removed: T2T-CHM13")
        self.assertNotRegex(" ".join(labels.values()), re.compile("phix", re.I))

    def test_no_file_claims_phix_removal(self):
        offenders = []
        for path in sorted(ROOT.rglob("*")):
            rel = path.relative_to(ROOT).as_posix()
            if (not path.is_file() or rel.startswith(".git/") or rel in EXPLAINERS
                    or path.suffix not in (".md", ".py", ".nf", ".config", ".html", ".R", ".r", ".sh")):
                continue
            if CLAIM.search(path.read_text(errors="replace")):
                offenders.append(rel)
        self.assertEqual(offenders, [])

    def test_reference_doc_and_links(self):
        doc = (ROOT / "docs" / "host_depletion_references.md").read_text()
        self.assertIn("there is no phiX in it", doc)
        self.assertIn("25 sequences", doc)
        for path in ("README.md", "docs/parameters.md", "docs/outputs.md"):
            with self.subTest(path=path):
                self.assertIn("host_depletion_references.md", (ROOT / path).read_text())

    def test_names_unchanged(self):
        """Documented, not changed: the parameter and the QC columns keep their names."""
        self.assertIn("params.t2t_phix_db", (ROOT / "main.nf").read_text())
        qc = (ROOT / "scripts" / "build_qc_summary.py").read_text()
        self.assertIn('"reads_after_t2t_phix"', qc)
        self.assertIn('"pct_removed_t2t_phix"', qc)

    def test_the_build_script_never_downloads_phix(self):
        """The code-level reason the index has no phiX. If someone fixes the download, this fails so
        the docs (and conf/databases.tsv, and the 'not removed' statements) are revisited."""
        script = (ROOT / "scripts" / "download_references.sh").read_text()
        code = "\n".join(l for l in script.splitlines() if not l.lstrip().startswith("#"))
        self.assertIn("ref/NC_001422.fna > ref/human-GCA-phix.fna", code)
        self.assertNotRegex(code, r"(wget|curl)[^\n]*NC_001422")
        self.assertNotIn("set -e", code)


if __name__ == "__main__":
    unittest.main()
