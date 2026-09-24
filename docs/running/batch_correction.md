# Batch correction (opt-in)

`--run_batch_correction true` removes batch effects from the decontaminated table with ConQuR
(`scripts/2500703_batch_correction_normalization.r`), then normalises it and measures how much
variance batch explains before and after. It takes its input from decontamination, or from
consensus/Bracken when decontamination is off, or from an existing table with
`--start_from batch_correction --decontam_otu_table <path> --metadata_file <path>`. It uses
`--batch_column`, `--sample_id_column` and `--type_column` from the decontamination settings.

| Parameter | Default | Meaning |
|---|---|---|
| `--batch_corr_covariates` | `age_diag,sex,bmi` | comma-separated metadata columns ConQuR adjusts for. Every named column must exist in the metadata, or the run stops with a named error |
| `--batch_corr_method` | `tune` | `tune` (auto-phased Tune_ConQuR) or `vanilla` (default ConQuR settings). Nothing else is accepted: ComBat is not implemented |
| `--phase` | `auto` | only `auto` is implemented: phase 1 (16 settings), then phase 3 if phase 1 did not reduce the batch effect enough |
| `--r2_threshold` | `0.25` | phase 1 is accepted if it cuts the batch R² by at least this fraction |
| `--tumor_only` | `false` | correct tumours only. Which types count as tumours is `--batch_corr_tumor_values` |
| `--batch_corr_tumor_values` | `Tumor` | comma-separated, case-insensitive. Not the same as decontamination's `--positive_values`, which means "not a control" |
| `--batch_corr_allow_uncorrected` | `false` | a failed correction stops the run. `true` continues with the raw counts, published as `corrected/UNCORRECTED_passthrough.tsv`, never under a corrected name |
| `--batch_corr_env` | pre-built TSCC env | see *Conda/Mamba* |
| `--conqur_sha` | `ff233085…` | the ConQuR commit the script was written against; any other install is refused |

**What ConQuR needs**, and what happens without it:

- **At least one covariate that varies within every batch.** `--batch_corr_covariates none`
  fails loudly, by design: ConQuR builds its model from covariates plus batch and cannot run on
  batch alone. A covariate that is constant within any one batch is dropped with a warning; if
  that leaves none, the tuned method fails.
- **At least 20 samples** for the tuned method, after these removals, each logged: samples that
  are all-zero after decontamination (listed in `batch_correction_sample_status.tsv`), samples
  with an unknown or missing batch, samples missing any covariate value, and batches with fewer
  than 5 samples.

A failure writes `corrected/BATCH_CORRECTION_STATUS.tsv` (`failed`, the reason and sample counts)
and stops the run, unless `--batch_corr_allow_uncorrected true`.

**Outputs**, under `05_BATCH_CORRECTION/<run>/`: `corrected/ConQuR_tuned.tsv` and
`ConQuR_tuned_parameters.txt` (the file keeps this name under `--batch_corr_method vanilla`
too; its parameters file then reads `Phase used: vanilla`), `corrected/BATCH_CORRECTION_STATUS.tsv`,
`normalized/` (the raw and corrected tables under several transforms), `pcoa_plots/`,
`permanova_summary.tsv` (batch and covariate R² for every table and transform; check this to
judge the correction) and `batch_correction_sample_status.tsv`.

---
[Back to the README](../../README.md)
