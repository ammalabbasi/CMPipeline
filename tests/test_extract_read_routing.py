#!/usr/bin/env python3
"""Guard the read-extraction routing in Modules/extract_reads.nf.

`samtools fastq -N -o FILE` writes only READ1/READ2 to FILE and sends category-0
records -- reads flagged neither READ1 nor READ2 -- to STDOUT. Inside a Nextflow task
stdout is .command.out, so that form both dropped those reads from the extracted FASTQ
and wrote read sequences into the task log.

Verified on samtools 1.21 against a synthetic BAM and CRAM carrying flags 4, 77, 141,
orphan mates and category-0 records: `-o FILE` wrote 3 of 5 reads and leaked the other
2 to stdout; naming neither -o nor -0 wrote all 5 and leaked nothing. Same defect as
pksProfiler F01.

These are static checks on the module text -- they need no samtools and no data. They
exist so the routing cannot silently regress.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "Modules" / "extract_reads.nf"


class ExtractReadRoutingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = MODULE.read_text()
        # the samtools fastq invocation, comments stripped
        cls.code = "\n".join(
            line for line in cls.text.splitlines() if not line.strip().startswith("#")
        )

    def test_no_output_routing_flags(self):
        """-o/-0/-1/-2/-s route categories to separate destinations; none may return."""
        for flag in ("-o ", "-0 ", "-1 ", "-2 ", "-s "):
            self.assertNotIn(
                f"samtools fastq" + flag, self.code.replace("\n", " "),
                f"samtools fastq must not use {flag.strip()}",
            )
        fastq_calls = [
            ln for ln in self.code.splitlines() if "samtools fastq" in ln
        ]
        self.assertTrue(fastq_calls, "no samtools fastq call found")
        for call in fastq_calls:
            for flag in (" -o ", " -0 ", " -1 ", " -2 ", " -s "):
                self.assertNotIn(
                    flag, call,
                    f"category routing flag {flag.strip()} reappeared in: {call.strip()}",
                )

    def test_never_discards_a_category(self):
        """No read category may be routed to /dev/null.

        Scoped to the extraction pipeline itself: `2>/dev/null` elsewhere is stderr
        suppression (htsfile probing the format) and is unrelated.
        """
        for line in self.code.splitlines():
            if "samtools fastq" in line or "samtools view" in line:
                self.assertNotIn(
                    "/dev/null", line,
                    f"a read category is being discarded: {line.strip()}",
                )

    def test_stream_is_compressed_not_left_on_stdout(self):
        """The fastq stream must be piped into a compressor, not written by samtools."""
        self.assertRegex(
            self.code.replace("\\\n", " "),
            r"samtools fastq[^|]*\|\s*\n?\s*bgzip",
            "samtools fastq output must be piped into bgzip",
        )

    def test_keeps_read_name_suffixes(self):
        self.assertIn("-N", self.code, "-N keeps /1 /2 suffixes on read names")

    def test_keeps_unmapped_selection_and_excludes_secondary_supplementary(self):
        self.assertIn("-f 4", self.code, "must select unmapped reads")
        self.assertIn("-F 2304", self.code, "must exclude secondary and supplementary")

    def test_keeps_cram_reference_arguments(self):
        self.assertIn(
            "REFERENCE_ARGS", self.code,
            "CRAM reference arguments must still reach samtools view",
        )

    def test_output_is_verified_readable(self):
        self.assertIn(
            "gzip -t", self.code,
            "the written FASTQ must still be integrity-checked",
        )


if __name__ == "__main__":
    unittest.main()
