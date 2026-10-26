# Outputs

Results are written under `--outdir`, which defaults to `RESULTS/` in the directory you launch
from, next to that run's `work/` and `pipeline_info/`:

| Directory | Contents |
|---|---|
| `UNMAPPED_BAM` | `<sample>.reads_after_filtering.txt`, `<sample>.extract_counts.tsv` and `<sample>.upstream_aligner.tsv` (BAM/CRAM) and the fastp JSON/HTML reports; the extracted and filtered FASTQ, and `<sample>.low_as_read_ids.txt.gz` (rescued read names), only with `--save_intermediates true` |
| `FASTQC` | per-stage FastQC reports (raw, filtered, post-GRCh38, post-T2T, post-pangenome) |
| `MAPPED_READS` | `<sample>.host_depleted.fastq.gz` and `<sample>.depletion_counts.tsv`; `.hg38.fastq.gz` and `.hg38.t2t.fastq.gz` only with `--save_intermediates true` |
| `BRACKEN` | KrakenUniq reports, Bracken genus/species tables, merged tables |
| `METAPHLAN4` | per-sample profiles and the merged abundance tables |
| `CONSENSUS_TAXA` | the consensus genus and species tables, `consensus_status.tsv`, `bracken_library_sizes.tsv`, and the Bracken vs MetaPhlAn prevalence plot |
| `HUMANN3` | gene families, pathway abundance and coverage, plus `merged/` |
| `04_DECONTAMINATION`, `05_BATCH_CORRECTION`, `06_VALIDATION` | optional stages; see their sections |
| `MULTIQC` | `multiqc_report.html` over fastp JSON and all five FastQC passes |
| `QC_SUMMARY` | `cmp.qc.summary.tsv` (one row per sample, 48 columns) and `cmp_cohort_report.html` |

Each directory can be moved on its own with the matching parameter: `--unmapped_bam_dir`,
`--fastqc_dir`, `--mapped_reads_dir`, `--krakenuniq_bracken_dir`, `--metaphlan4_dir`,
`--consensus_taxa_dir`, `--humann3_dir`, `--decontam_dir`, `--batch_corr_dir`, `--validation_dir`,
`--multiqc_dir` and `--qc_summary_dir`. Each defaults to its name under `--outdir`.

Run records are written **outside** the repository, under the launch directory, so
concurrent runs do not overwrite each other's measurements:

| `pipeline_info/` | Contents |
|---|---|
| `RUN_REPORT.txt` | the run in reading order, appended per session — status, reproducibility, per-process cost, the conda envs specified, and *TOOL VERSIONS (as run)*: what each task actually executed |
| `provenance/<sessionId>_<runName>.json` | the complete machine-readable record, one per launch (a `-resume` adds a new file). It carries SHA-256 digests of the cohort tables (`outputs.cohort_tables`), the conda envs this run used, ConQuR/cqrReg versions, and a digest of `conda_envs/` |
| `trace.txt`, `report.html`, `timeline.html`, `dag.html` | Nextflow's own execution records, on by default ; `trace.txt` adds `attempt`, `cpus`, `memory` and `workdir` to Nextflow's default columns |
| `<work dir>/.tool_versions.tsv` | per task: each tool the task ran, its package, status, version, lockfile pin, resolved path and the active `$CONDA_PREFIX` (`scripts/record_tool_versions.sh`) |

**What ran, not just what was pinned.** The *CONDA ENVS THIS RUN USED* section lists what each
lockfile specifies. *TOOL VERSIONS (as run)* lists what each task's binaries reported from inside
the task, after conda activation: tool, version, pin, the version conda installed in the prefix,
task count and where it ran from. Anything that is not the pinned build from the task's own env
appears under *TOOL VERSION MISMATCH*: `OUTSIDE_ENV` (resolved outside the env, e.g. a login
`tools/minimap2` earlier on PATH), `VERSION_MISMATCH`, `NOT_FOUND`, `NO_ENV`, `MISSING` (a recording
task left no record) and `INSTALLED_MISMATCH` (the binary disagrees with its env's conda-meta). With
`--tool_version_check strict` (default) a task stops the run at the first of the first four, and the
finaliser exits non-zero on any of them; the JSON record carries the same under
`environment.tool_versions`. Tool names, versions and paths only, never a sample ID.

`RUN_REPORT.txt` states the commit and **whether the working tree was dirty**. A run
launched with uncommitted edits says so, in those words — "this run is NOT attributable to
that commit" — because it is not.

`cmp.qc.summary.tsv` is the one place that answers "what happened to this sample": reads in
and out of filtering, retention at each depletion stage, whether it cleared
`--consensus_min_reads`, how many genera and species each profiler called and how many they
agreed on, its decontamination role, and whether batch correction included it at each rank
(`batch_species_status`, `batch_genus_status`: `included`, or `excluded: <reason>`). A run that
starts at `--start_from decontam` or `batch_correction` still gets one row per analysed sample,
from those stages; the upstream columns then read `NA`. `cleared_consensus_gate` reuses the exact
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

**Low-score mapped rescue (BAM/CRAM, `--extract_low_score_mapped`).** Ten columns record, per
sample, what the upstream aligner was and whether reads it reported as mapped with a low score were
also extracted (see [running/classification.md](running/classification.md#low-score-mapped-reads-tcga-gdc-bams)).
All `NA` for FASTQ input.

| Column | Meaning |
|---|---|
| `upstream_aligner`, `upstream_aligner_version` | the aligner's `@PG` PN and VN (bwa mem's own record when present; post-processing tools such as samtools or Picard are never reported); the command line is never recorded |
| `upstream_bwa_min_score` | the `-T` in the `bwa mem` command line: a number, `default` (not given, i.e. 30), `unparsed`; `NA` without bwa mem |
| `low_as_mode` | `on` or `off` for this sample |
| `low_as_mode_reason` | why: `auto_bwa_mem_T_below_30` (on), `auto_bwa_mem_T_default`, `auto_bwa_mem_T_30_or_more`, `auto_bwa_mem_T_unparsed`, `auto_bwa_not_mem`, `auto_other_aligner`, `auto_no_pg`, `auto_no_aligner_pg` (only samtools/Picard-style post-processing @PG), `auto_unknown_aligner` (all off), `requested_on`, `requested_off` |
| `low_as_reads_extracted` | primary mapped records with AS < 30 extracted in addition to the FLAG 4 reads (`0` when off) |
| `low_as_no_as_records` | primary mapped records without an AS tag (never extracted; counted only when the rescue is on, `NA` otherwise) |
| `low_as_reads_after_depletion` | rescued reads still present in the host-depleted FASTQ |
| `low_as_reads_classified` | rescued reads KrakenUniq classified |
| `low_as_share_of_bracken` | `low_as_reads_classified` as a percentage of `bracken_microbial_reads` (KrakenUniq classified reads, the input Bracken re-estimates from). Bracken does not attribute individual reads, so this is the rescued share of its input, not of each taxon |

With the rescue on, `extracted_unmapped_records` keeps its meaning (FLAG 4 only) and
`reads_raw` (the reads entering fastp) is `extracted_unmapped_records + low_as_reads_extracted`.
The per-sample survival rows are also written to `QC_SUMMARY/low_as_survival/<sample>.low_as_survival.tsv`,
and `RUN_REPORT.txt` gets a counts-only *LOW-SCORE RESCUE* section (inputs per mode and reason, upstream
aligners, rescued reads summed per stage).

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
- **`RESULTS/QC_SUMMARY/cmp.qc.summary.tsv`** — one row per sample, 48 columns (see *Outputs*).
- **`RESULTS/QC_SUMMARY/cmp_cohort_report.html`** — a self-contained page (no network, no
  CDN) with the retention chart and the tables. Disable with `--skip_cohort_report true`.
- **`pipeline_info/RUN_REPORT.txt`** and the provenance record — see *Outputs*.

The cohort report reads the summary table, not the results tree, so the two cannot
disagree about what a run produced.

---
[Back to the README](../README.md)
