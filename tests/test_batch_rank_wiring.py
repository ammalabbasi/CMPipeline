#!/usr/bin/env python3
"""Batch correction runs on both ranks wherever a species-first route produces both.

Static checks on main.nf and the module; the executed wiring is checked by `nextflow run
-preview` and by the CRC50 downstream run.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class BatchRankWiring(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.main = (ROOT / "main.nf").read_text()
        cls.module = (ROOT / "Modules" / "batch_correction.nf").read_text()

    def test_every_route_tags_the_rank(self):
        for fragment in (
            "tuple(prefix, table, metadata, 'species')",                  # decontam species
            "aggregateDecontaminatedSpecies.out.genus",                   # decontam genus roll-up
            "tuple('validated', t, validatedMeta, 'genus')",              # validation, no decontam
            "tuple('consensus_run', t, consensusMeta, 'genus')",          # consensus, no decontam
            "tuple('bracken_run', bracken_genus_file, brackenMeta, 'genus')",  # Bracken, no decontam
            "params.taxonomic_rank\n            )).set { ch_batch_input }",  # supplied table
        ):
            self.assertIn(fragment, self.main)

    def test_no_route_selects_a_single_rank_any_more(self):
        self.assertNotIn("CONSENSUS_OUTPUT", self.main)
        self.assertNotIn("params.taxonomic_rank == 'species' ? VALIDATED_SPECIES", self.main)

    def test_module_keeps_ranks_apart(self):
        self.assertIn("tuple val(prefix), path(otu_table), path(metadata), val(rank)", self.module)
        self.assertIn('"${params.batch_corr_dir}/${prefix}_${rank}"', self.module)
        self.assertIn('--prefix "${prefix}_${rank}"', self.module)
        self.assertIn('path("${rank}.batch_sample_status.tsv"), emit: qc_status', self.module)

    def test_batch_ledger_reaches_the_qc_summary(self):
        self.assertIn("qcParts << BatchCorrection.out.qc_status", self.main)


if __name__ == "__main__":
    unittest.main()
