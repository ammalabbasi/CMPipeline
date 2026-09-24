# Functional profiling (HUMAnN 3)

HUMAnN3 is off by default (`--skip_humann3 true`). When enabled it needs
`--humann3_nucleotide_db` and `--humann3_protein_db`, both checked before the run starts:

```bash
nextflow run main.nf --skip_humann3 false \
  --humann3_nucleotide_db /db/chocophlan --humann3_protein_db /db/uniref
```

**HUMAnN3 runs only on samples that cleared `--consensus_min_reads`**, the same samples
MetaPhlAn ran on, and each is handed its own MetaPhlAn profile, so HUMAnN3 does not recompute
one. Samples below the gate are skipped for functional profiling: HUMAnN3 would find almost no
species for its nucleotide search there and fall back to a translated search of every read, its
slowest mode and the one most exposed to residual human reads. Which samples ran is visible in
the QC summary (`metaphlan_genera` is `NA` for a skipped sample). If no sample cleared the gate,
the run logs a warning and produces no functional profiles. The profile is passed as a
process input rather than read from `RESULTS/METAPHLAN4/`, so the result does not depend
on which task finished first.

---
[Back to the README](../../README.md)
