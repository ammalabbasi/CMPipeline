# Independent batch environment and reproducible execution

These are package-supported commands, not user-side R patches. The supplied batch
lockfile targets Linux x86_64. Run installation and analysis on a compute node when
using HPC; obtain sufficient resources with the site's normal scheduler commands.

## Fresh user-owned environment

From the repository root:

```bash
bash scripts/setup_batch_env.sh --prefix /absolute/new/batch-env --conda /path/to/conda
```

The prefix must not exist. The helper installs the exact batch lockfile, then
checksum-verified cqrReg 1.2.1 and ConQuR_par commit
`ff233085cedc36a24382038067a4cece3e94e0a0`. It does not overwrite a shared prefix or
install packages during analysis. A failed partial prefix is retained for diagnosis;
use a new prefix for an explicitly approved retry. No automatic deletion/rebuild.

For an offline R-source cache, supply `--source-dir /path/to/tarballs`, containing
`cqrReg_1.2.1.tar.gz` and
`ConQuR_par-ff233085cedc36a24382038067a4cece3e94e0a0.tar.gz`. The helper checks its
embedded expected hashes, not a cache owner's mutable manifest. Conda packages must
still be cached or downloadable. Without `--source-dir`, the helper downloads these
sources from CRAN and GitHub and checks the same hashes.

`cmp_batch_versions.tsv`, `cmp_batch_lock.sha256`, and `cmp_batch_sources.sha256`
inside the new prefix document setup; the normal pipeline provenance records the
environment actually used. Source-download and package-cache availability are
external prerequisites; an offline-cache test does not certify online retrieval.

## Official standalone batch entry

```bash
nextflow run main.nf -profile local -c my-resources.config \
  --start_from batch_correction --run_batch_correction true --run_decontam false \
  --decontam_otu_table counts.tsv --metadata_file metadata.csv \
  --sample_id_column patient --type_column sample_type --batch_column shipment_batch \
  --taxon_column clade_name --taxonomic_rank species \
  --batch_corr_covariates sample_type --batch_corr_seed 314159 \
  --batch_corr_env /absolute/new/batch-env --outdir /absolute/new/results
```

Choose covariates from the actual study design; `sample_type` here illustrates a
functional Tumor/Normal test, not a replacement for required age/sex/BMI adjustment.
The tuned minimum of 20 eligible samples and per-batch selection rules are unchanged.
The portable default no longer assumes a TSCC batch prefix. The `tscc` profile keeps
the site's optional shared prefix; an explicit `--batch_corr_env` overrides it.

For a bounded two-CPU test allocation, `my-resources.config` can contain:

```groovy
executor.cpus = 2
executor.memory = '16.GB'
process {
    withName: BatchCorrection { cpus = 2; memory = '8.GB'; time = '45m' }
}
```

Set production resources appropriately; this is not a recommendation for arbitrary
cohort sizes. Batch workers use `task.cpus`; BLAS/OpenMP threads are capped at one.
Use the same CPU count, seed (default 1), input ordering, code and environment for
independent comparisons. R uses L'Ecuyer-CMRG/Inversion/Rejection and Unix multicore
worker streams. Cross-CPU-count or cross-platform identity is not promised. R runs
with `--vanilla` and environment-local libraries to avoid personal startup/library
contamination. The original ConQuR dependency/parameter search is not rewritten.

## Reporting: execution is not efficacy

- `batch_correction_run_settings.tsv`: seed, worker count, RNG kinds and model fields.
- `permanova_summary.tsv`: evaluation mode, distance, model and N for each result.
- `batch_correction_quality.tsv`: observed R2 values, raw diagnostic relative change,
  interpretable percentage (NA for nonpositive baseline/negative final/nonfinite
  metrics), execution status and explicit unassessed biological-signal preservation.
- `batch_correction_methods.tsv`: selected original/correction method per prevalence
  bin and count of taxa represented. Vanilla/explicit uncorrected modes are labelled.

An unchanged result may be an intentional original-data selection, a negligible
adjustment, or an unobserved per-taxon fit problem. The pinned dependency does not
expose every taxon's fit outcome; reports say `not_exposed_by_pinned_dependency`
rather than inventing one. Negative R2 is preserved, not clamped. Same-cohort tuning
statistics are descriptive and do not certify biological-signal preservation.
Internal phase-selection rules are unchanged; noninterpretable reporting does not
silently choose a different correction. Packaging/scientific limitations must retain
these boundaries until separately addressed and tested.
