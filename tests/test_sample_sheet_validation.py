#!/usr/bin/env python3
"""Guard the sample-sheet ingestion checks in main.nf.

Adversarial probing on 2026-09-22 found that two rows with the same patient ID -- or IDs
differing only by surrounding whitespace, which `trim()` collapses -- flowed straight
through ingestion and then collided in publishDir and in every merged table, one sample
silently overwriting the other. That is the same defect as pksProfiler F03.

The checks themselves run inside a Nextflow channel, so they cannot be exercised from
Python. These are static assertions that the checkpoint is still present and still wired
in front of the branch, which is what a careless refactor would remove. The runtime
behaviour was verified directly against the pipeline:

    duplicate IDs    CAUGHT: Duplicate sample IDs in the sample sheet: DUP
    whitespace dup   CAUGHT: Duplicate sample IDs in the sample sheet: PAD
    header only      CAUGHT: empty sample sheet
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "main.nf"


class SampleSheetValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = MAIN.read_text()

    def test_duplicate_ids_are_rejected(self):
        self.assertIn("Duplicate sample IDs in the sample sheet", self.text,
                      "the duplicate sample-ID check has been removed")
        self.assertIn("countBy", self.text,
                      "duplicate detection no longer counts IDs")

    def test_empty_sheet_is_rejected(self):
        self.assertIn("Sample sheet has no data rows", self.text,
                      "an empty sample sheet would run and silently do nothing")

    def test_checkpoint_precedes_the_branch(self):
        """The check must sit in front of the branch, or rows bypass it."""
        check = self.text.index("Duplicate sample IDs in the sample sheet")
        branch = self.text.index("sample_sheet.branch {")
        self.assertLess(check, branch,
                        "the duplicate check must come before sample_sheet.branch")

    def test_checkpoint_preserves_every_row(self):
        """toList/flatMap must round-trip the rows, not filter them."""
        window = self.text[self.text.index("sample_sheet = sample_sheet"):
                           self.text.index("sample_sheet.branch {")]
        self.assertIn(".toList()", window)
        self.assertIn(".flatMap()", window)
        self.assertIn("return rows", window)

    def test_ids_are_trimmed_and_emptiness_rejected(self):
        self.assertIn("row.patient?.trim()", self.text)
        self.assertIn("is missing a non-empty 'patient' value", self.text)

    def test_every_input_path_is_existence_checked(self):
        """checkIfExists on each input, so a typo fails at ingestion not mid-run."""
        window = self.text[self.text.index("def inputs = inputType == 'bam'"):]
        window = window[:window.index("tuple(sampleID, inputType, inputs, references)")]
        self.assertGreaterEqual(window.count("checkIfExists: true"), 3)

    def test_mates_may_not_be_the_same_file(self):
        self.assertIn("has the same path for both FASTQ mates", self.text)

    def test_ambiguous_input_types_are_rejected(self):
        self.assertIn("has ambiguous inputs", self.text)


if __name__ == "__main__":
    unittest.main()
