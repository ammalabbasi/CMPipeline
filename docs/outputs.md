# Outputs

Results are written under `--outdir`, which defaults to `RESULTS/` in the directory you launch
from, next to that run's `work/` and `pipeline_info/`:

| Directory | Contents |
|---|---|
| `UNMAPPED_BAM` | `<sample>.reads_after_filtering.txt`, `<sample>.extract_counts.tsv` (BAM/CRAM) and the fastp JSON/HTML reports; the extracted and filtered FASTQ only with `--save_intermediates true` |
| `FASTQC` | per-stage FastQC reports (raw, filtered, post-GRCh38, post-T2T, post-pangenome) |
| `MAPPED_READS` | `<sample>.host_depleted.fastq.gz` and `<sample>.depletion_counts.tsv`; `.hg38.fastq.gz` and `.hg38.t2t.fastq.gz` only with `--save_intermediates true` |
| `BRACKEN` | KrakenUniq reports, Bracken genus/species tables, merged tables |
| `METAPHLAN4` | per-sample profiles and the merged abundance tables |
| `CONSENSUS_TAXA` | the consensus genus and species tables, `consensus_status.tsv`, `bracken_library_sizes.tsv`, and the Bracken vs MetaPhlAn prevalence plot |
| `HUMANN3` | gene families, pathway abundance and coverage, plus `merged/` |
| `04_DECONTAMINATION`, `05_BATCH_CORRECTION`, `06_VALIDATION` | optional stages; see their sections |
| `MULTIQC` | `multiqc_report.html` over fastp JSON and all five FastQC passes |
| `QC_SUMMARY` | `cmp.qc.summary.tsv` (one row per sample, 36 columns) and `cmp_cohort_report.html` |

Each directory can be moved on its own with the matching parameter: `--unmapped_bam_dir`,
`--fastqc_dir`, `--mapped_reads_dir`, `--krakenuniq_bracken_dir`, `--metaphlan4_dir`,
`--consensus_taxa_dir`, `--humann3_dir`, `--decontam_dir`, `--batch_corr_dir`, `--validation_dir`,
`--multiqc_dir` and `--qc_summary_dir`. Each defaults to its name under `--outdir`.

Run records are written **outside** the repository, under the launch directory, so
concurrent runs do not overwrite each other's measurements:

| `pipeline_info/` | Contents |
|---|---|
| `RUN_REPORT.txt` | the run in reading order, appended per session — status, reproducibility, per-process cost |
| `provenance/<sessionId>_<runName>.json` | the complete machine-readable record, one per launch (a `-resume` adds a new file). It carries SHA-256 digests of the cohort tables (`outputs.cohort_tables`), the conda envs this run used, ConQuR/cqrReg versions, and a digest of `conda_envs/` |
| `trace.txt`, `report.html`, `timeline.html`, `dag.html` | Nextflow's own execution records, on by default |

`RUN_REPORT.txt` states the commit and **whether the working tree was dirty**. A run
launched with uncommitted edits says so, in those words — "this run is NOT attributable to
that commit" — because it is not.

`cmp.qc.summary.tsv` is the one place that answers "what happened to this sample": reads in
and out of filtering, retention at each depletion stage, whether it cleared
`--consensus_min_reads`, how many genera and species each profiler called and how many they
agreed on, and its decontamination role. `cleared_consensus_gate` reuses the exact
definition `Bracken.nf` gates on, so the table can never disagree with the pipeline.

**Read retention, and what it is relative to.** For FASTQ input, `reads_raw` is the whole
library. For BAM/CRAM input the pipeline extracts only the reads the aligner left unmapped, so
`reads_raw` is that extracted subset, not the library. Three columns carry the library for
BAM/CRAM (all `NA` for FASTQ):

| Column | Meaning |
|---|---|
| `library_primary_records` | primary records in the BAM/CRAM, the library size |
| `extracted_unmapped_records` | how many of them were unmapped and extracted; `pct_library_unmapped` is this as a percentage of the library |
| `pct_retained_of_library` | host-depleted reads as a percentage of the whole library |

`pct_retained_overall` is host-depleted reads as a percentage of the reads **entering fastp**:
the library for FASTQ, but only the extracted unmapped subset for BAM/CRAM. Compare samples of
different input types with `pct_retained_of_library`, not `pct_retained_overall`.

`sample_status` says in one word what happened to each sample: `dropped_zero_after_fastp`,
`no_reads_after_depletion`, `below_bracken_threshold` (Bracken made no call, and `bracken_status`
is `no_call`) or `profiled`. A sample in the sheet that produced nothing (a failed task, or one dropped by `--sample_failure_strategy ignore`) still gets a row, with `sample_status` `no_outputs`.

`<sample>.reads_after_filtering.txt` holds the post-fastp read count; it is what the
workflow uses to decide whether a sample continues, so it is also the first place to look
when a sample is missing from downstream output.

Work files remain in `work/` and should be retained when using `-resume`. Final tables and
reports are copied out of `work/` into `--outdir`; intermediate read files are copied only with
`--save_intermediates true`, so by default they exist in `work/` alone.

## Quality control and reporting

Every run produces, unless switched off:

- **`RESULTS/MULTIQC/multiqc_report.html`** — fastp JSON plus all five FastQC passes
  (raw, filtered, post-GRCh38, post-T2T, post-pangenome). Disable with `--skip_multiqc true`.
- **`RESULTS/QC_SUMMARY/cmp.qc.summary.tsv`** — one row per sample, 36 columns (see *Outputs*).
- **`RESULTS/QC_SUMMARY/cmp_cohort_report.html`** — a self-contained page (no network, no
  CDN) with the retention chart and the tables. Disable with `--skip_cohort_report true`.
- **`pipeline_info/RUN_REPORT.txt`** and the provenance record — see *Outputs*.

The cohort report reads the summary table, not the results tree, so the two cannot
disagree about what a run produced.

---
[Back to the README](../README.md)
