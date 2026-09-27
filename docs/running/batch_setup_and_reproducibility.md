# Batch environment and reproducible execution

The supplied batch lockfile targets Linux x86_64. Run analysis on a compute node when using HPC;
obtain sufficient resources with the site's normal scheduler commands.

## Environment

Nextflow builds the batch-correction environment from
`conda_envs/lock/batch_correction_env.linux-64.txt`, exactly like the other environments; nothing
needs to be installed beforehand. The lockfile ends with two packages that no conda channel
provides, r-cqrreg 1.2.1 and r-conqur (`ivartb/ConQuR_par` at
`ff233085cedc36a24382038067a4cece3e94e0a0`). They are committed as conda packages in
`conda_envs/local_channel/` and named in the lockfile under `@LOCAL_CHANNEL@`, which the workflow
replaces with that directory's absolute path (the resolved lockfile is written to
`<work-dir>/conda_locks/`, under a name that includes a digest of its content, and is never
rewritten; an existing file whose content differs stops the run). The file is created by an atomic
move; on a file system without atomic moves that step fails with an error rather than continuing. The path goes into the lockfile unescaped, so the pipeline checkout path may contain
only letters, digits and `. _ / + -`; anything else (including spaces, `@` and a literal `~`)
stops the run with that message.
Moving the checkout changes the lockfile, so Nextflow builds the environment again.
The conda cache defaults to `<checkout>/.conda_cache` (`nextflow.config`): for a read-only
checkout set `NXF_CONDA_CACHEDIR` to a writable directory. They were built by `scripts/build_local_conda_packages.sh` from
`conda_envs/recipes/`, from the same checksum-verified sources and against the lockfile's own
r-base and compilers; the lockfile pins them by md5. The upstream sources and license texts are
committed too (`conda_envs/local_channel/NOTICE.md`). Before running, the module checks that the
installed ConQuR is the pinned commit.

Optional: `bash scripts/setup_batch_env.sh --prefix /absolute/new/batch-env` builds the same
environment into a prefix for `--batch_corr_env` (e.g. to reuse it across runs). The prefix must
not exist; `cmp_batch_versions.tsv` and `cmp_batch_lock.sha256` inside it record the build.

## Official standalone batch entry

```bash
nextflow run main.nf -profile local -c my-resources.config \
  --start_from batch_correction --run_batch_correction true --run_decontam false \
  --decontam_otu_table counts.tsv --metadata_file metadata.csv \
  --sample_id_column patient --type_column sample_type --batch_column shipment_batch \
  --taxon_column clade_name --taxonomic_rank species \
  --batch_corr_covariates sample_type --batch_corr_seed 314159 \
  --outdir /absolute/new/results
```

Choose covariates from the actual study design; `sample_type` here illustrates a
functional Tumor/Normal test, not a replacement for required age/sex/BMI adjustment.
The tuned minimum of 20 eligible samples and per-batch selection rules are unchanged.
No profile needs a site-specific batch prefix; `--batch_corr_env` is only an optional override.

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
