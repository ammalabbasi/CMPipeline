# All parameters

Every `--parameter` `main.nf` defines, grouped as in `main.nf`. Booleans accept only `true`/`false`.
Paths and values are validated before the first task runs. The run-mode guides in
[`docs/running/`](running/) explain each group in context.

## Workflow Control Parameters

| Parameter | Default | Description |
|---|---|---|
| `--skip_host_depletion` | `false` | skip the GRCh38/T2T/pangenome depletion stage |
| `--skip_classification` | `false` | skip KrakenUniq/Bracken and MetaPhlAn |
| `--skip_consensus` | `false` | skip consensus (Bracken ∩ MetaPhlAn agreement); Bracken tables then go downstream |
| `--skip_humann3` | `true` | HUMAnN 3 is off by default; `false` turns it on (needs both HUMAnN databases) |
| `--run_decontam` | `true` | run decontamination (on by default; needs `--metadata_file` and the four role parameters) |
| `--run_batch_correction` | `false` | run ConQuR batch correction (off by default; needs covariates and >= 20 samples) |
| `--skip_multiqc` | `false` | skip the MultiQC report (fastp JSON + all FastQC passes) |
| `--skip_cohort_report` | `false` | skip the self-contained HTML cohort report |
| `--sample_failure_strategy` | `finish` | what one sample's failure does. `finish` stops the run once running tasks end; `ignore` lets the other samples complete and lists the failures in `RUN_REPORT.txt` and `trace.txt`. Cohort-level steps (merges, consensus, decontamination, batch correction, reports) always stop on failure |
| `--save_intermediates` | `false` | also publish intermediate read files (extracted, filtered and per-stage depleted FASTQ, MetaPhlAn Bowtie2 output). By default only the host-depleted FASTQ is published; everything else stays in `work/` |
| `--start_from` | `beginning` | `decontam` or `batch_correction` enters at that stage from an existing table (see those sections) |
| `--consensus_otu_table` | `null` | the consensus table to decontaminate, with `--start_from decontam` |
| `--decontam_otu_table` | `null` | the decontaminated table to batch-correct, with `--start_from batch_correction` |
| `--metadata_file` | none (required) | per-sample metadata, tab- or comma-separated |

## Decontamination Parameters

| Parameter | Default | Description |
|---|---|---|
| `--decontam_threshold` | `0.1` | decontam's prevalence-test threshold |
| `--decontam_min_prevalence` | `0.05` | after that filter, a taxon must be present in at least 5% of samples to be tested |
| `--decontam_min_abundance` | `5` | a taxon's reads in a sample below this count as absent |
| `--decontam_min_batches` | `2` | a taxon flagged in at least this many batches is a contaminant |
| `--decontam_zero_total_policy` | `keep` | `keep`: an all-zero sample stays in the prevalence test, as long as its library was not empty. `drop`: every all-zero sample is excluded (the older behaviour, which could remove a batch's controls) |
| `--decontam_library_sizes` | none | the library-size table for `--start_from decontam`: pass the `bracken_library_sizes.tsv` that consensus wrote. In a full run it is supplied automatically |
| `--decontam_unevaluable_policy` | **`skip`** | `skip`: record the taxon as `unevaluable` in that batch and not flagged there; the `--decontam_min_batches` rule then decides from the batches that could test it. `error`: stop the run |
| `--decontam_min_library` | `1000` | under `keep`, an all-zero sample whose total Bracken library is below this is excluded as `insufficient_library` |
| `--decontam_allow_unknown_library` | `false` | an all-zero sample whose library size is unknown stops the run unless this is `true`, in which case it is kept |
| `--min_total_per_batch` | `5`, `1`, `1` | a batch with fewer samples, positives or controls than this is not tested |
| `--min_positive_per_batch` | `5`, `1`, `1` | a batch with fewer samples, positives or controls than this is not tested |
| `--min_control_per_batch` | `5`, `1`, `1` | a batch with fewer samples, positives or controls than this is not tested |
| `--batch_column` | `shipment_batch` | metadata column holding the batch; also used by batch correction |
| `--batch_corr_tumor_values` | `Tumor` | comma-separated, case-insensitive. Not the same as decontamination's `--positive_values`, which means "not a control" |
| `--batch_corr_allow_uncorrected` | `false` | a failed correction stops the run. `true` continues with the raw counts, published as `corrected/UNCORRECTED_passthrough.tsv`, never under a corrected name |
| `--conqur_sha` | `ff233085…` | the ConQuR commit the script was written against; any other install is refused |
| `--taxon_column` | `clade_name` | taxonomy column in the OTU table; every other column is a sample |
| `--taxonomic_rank` | `null` | **required when decontamination or batch correction runs**: `genus` or `species`, the table to decontaminate/correct |
| `--sample_id_column` | `sampleid` | metadata column holding the sample IDs |
| `--unselected_type_policy` | `error` | `error` (default) stops on a sample whose type is in neither role list; `exclude` drops it |
| `--cohort_preset` | `null` | `tcga` sets the TCGA role contract (type column `Sample Type_x`, positives Primary Tumor + Solid Tissue Normal, controls Blood Derived Normal) |
| `--type_column` | `params.cohort_preset == 'tcga' ? 'Sample Type_x' : 'sample_t` | metadata column holding each sample's type/role |
| `--control_role` | `params.cohort_preset == 'tcga' ? 'surrogate_biological' : nu` | **required for decontamination**: `experimental_blank` (true blanks) or `surrogate_biological` (a sample class standing in for blanks) |
| `--positive_values` | `params.cohort_preset == 'tcga' ? 'Primary Tumor,Solid Tissue` | **required for decontamination**: comma-separated sample types that are real material (not controls) |
| `--control_values` | `params.cohort_preset == 'tcga' ? 'Blood Derived Normal' : nu` | **required for decontamination**: comma-separated sample types used as negative controls |

## Batch Correction Parameters

| Parameter | Default | Description |
|---|---|---|
| `--batch_corr_covariates` | `age_diag,sex,bmi` | comma-separated metadata columns ConQuR adjusts for. Every named column must exist in the metadata, or the run stops with a named error |
| `--tumor_only` | `false` | correct tumours only. Which types count as tumours is `--batch_corr_tumor_values` |
| `--phase` | `auto` | only `auto` is implemented: phase 1 (16 settings), then phase 3 if phase 1 did not reduce the batch effect enough |
| `--r2_threshold` | `0.25` | phase 1 is accepted if it cuts the batch R² by at least this fraction |
| `--bracken_read_length` | `50` | Bracken `-r`; the Kraken database must hold a `database<N>mers.kmer_distrib` for it, checked at launch |
| `--bracken_threshold` | `2` | Bracken `-t`: the minimum reads a taxon needs to be re-estimated. A read threshold, not a thread count. A thin sample in which no taxon reaches it gets an empty no-call table instead of failing the run |
| `--consensus_min_reads` | `100000` | the MetaPhlAn gate above; a non-negative integer (`1e5` is rejected) |

## Input/Output Paths

| Parameter | Default | Description |
|---|---|---|
| `--sample` | `samples.csv` | **required**: the sample sheet CSV (`patient` + FASTQ, BAM or CRAM columns) |
| `--input_data_type` | `auto` | auto | bam | cram | fastq |
| `--cram_reference` | `null` | optional matching FASTA for CRAM decoding |
| `--outdir` | `RESULTS/` in the launch directory | where every stage publishes (see *Outputs*) |
| `--unmapped_bam_dir` | ``--outdir`/UNMAPPED_BAM` | where `UNMAPPED_BAM/` is published; defaults to `UNMAPPED_BAM` under `--outdir` |
| `--mapped_reads_dir` | ``--outdir`/MAPPED_READS` | where `MAPPED_READS/` is published; defaults to `MAPPED_READS` under `--outdir` |
| `--fastqc_dir` | ``--outdir`/FASTQC` | where `FASTQC/` is published; defaults to `FASTQC` under `--outdir` |
| `--multiqc_dir` | ``--outdir`/MULTIQC` | where `MULTIQC/` is published; defaults to `MULTIQC` under `--outdir` |
| `--qc_summary_dir` | ``--outdir`/QC_SUMMARY` | where `QC_SUMMARY/` is published; defaults to `QC_SUMMARY` under `--outdir` |
| `--conda_cache_dir` | `$NXF_CONDA_CACHEDIR`, else `.conda_cache/` in the repo | recorded in the run record; mirrors `conda.cacheDir` in `nextflow.config`, which is what Nextflow actually uses |
| `--krakenuniq_bracken_dir` | ``--outdir`/BRACKEN` | where `BRACKEN/` is published; defaults to `BRACKEN` under `--outdir` |
| `--metaphlan4_dir` | ``--outdir`/METAPHLAN4` | where `METAPHLAN4/` is published; defaults to `METAPHLAN4` under `--outdir` |
| `--humann3_dir` | ``--outdir`/HUMANN3` | where `HUMANN3/` is published; defaults to `HUMANN3` under `--outdir` |
| `--consensus_taxa_dir` | ``--outdir`/CONSENSUS_TAXA` | where `CONSENSUS_TAXA/` is published; defaults to `CONSENSUS_TAXA` under `--outdir` |
| `--decontam_dir` | ``--outdir`/04_DECONTAMINATION` | where `04_DECONTAMINATION/` is published; defaults to `04_DECONTAMINATION` under `--outdir` |
| `--batch_corr_dir` | ``--outdir`/05_BATCH_CORRECTION` | where `05_BATCH_CORRECTION/` is published; defaults to `05_BATCH_CORRECTION` under `--outdir` |
| `--antismash_dir` | `conda_envs/antismash_env.yml`, `ANTISMASH/` | reserved for the antiSMASH module in `Modules/antismash.nf`, which `main.nf` does not currently call |
| `--validation_dir` | ``--outdir`/06_VALIDATION` | where `06_VALIDATION/` is published; defaults to `06_VALIDATION` under `--outdir` |
| `--run_alignment_validation` | `false` | turn the lane on |
| `--validation_min_support` | `25` | reads on the best genome for a species to be listed, and for a sample x species to be scored at all |
| `--validation_min_confirmed` | 10 | a species with at least `--validation_min_support` Bracken reads but fewer supporting alignments than this is `unsupported`, and its count is zeroed (decision 29; PRISM's cutoff) |
| `--validation_min_bins_ratio` | `0.5` | 10 kb windows hit / windows expected for that many read pairs |
| `--validation_min_identity` | `0.98` | median alignment identity |
| `--validation_min_unique_frac` | `0.02` | reads unique to the species / all its reads |
| `--validation_max_genomes` | `50` | RefSeq genomes per species, before de-duplication |
| `--validation_dedup_ani` | `99.5` | skani ANI at which two genomes of a species count as one |
| `--validation_max_missing_frac` | `0.10` | the build fails if more than this fraction of species get no genome |
| `--validation_max_ref_gb` | 40 | refuse a planned reference larger than this many Gb of sequence, checked from RefSeq genome sizes **before** anything is downloaded or indexed |
| `--validation_batch_size` | `20` | samples per alignment task |
| `--validation_groups` | `ref/validation_groups.tsv` | species scored together as one taxon (E. coli + Shigella) |
| `--validation_export_reads` | `false` | also write each validated taxon's supporting reads (protected data) |
| `--validation_taxdb` | `<kraken_db>/taxDB` | the NCBI taxonomy used to match names to taxids |
| `--validation_cache_dir` | `.validation_cache/` in the repo, or `$CMP_VALIDATION_CACHE` | where built references are kept and reused |
| `--validation_ref` | none | a pre-built `ref-<key>/` folder from the cache: skips the species list and the build (no internet needed; pins one reference across cohorts) |
| `--validation_refseq_snapshot` | the month the run starts | RefSeq snapshot month (`YYYY-MM`), part of the reference's cache key. Within a month the cached reference is reused; a new month builds a fresh, correctly labelled one. A past month is reused only if still cached; otherwise the build stops (pin it with `--validation_ref`) |

## Databases And Reference Files

| Parameter | Default | Description |
|---|---|---|
| `--hg38_db` | `null` | Set in a site profile or with --hg38_db |
| `--t2t_phix_db` | `null` | Set in a site profile or with --t2t_phix_db |
| `--pangenome_db` | `null` | optional .mmi file or directory containing .mmi files |
| `--kraken_db` | `null` | Set in a site profile or with --kraken_db |
| `--metaphlan_db` | `null` | Set in a site profile or with --metaphlan_db |
| `--metaphlan_index` | `mpa_vJun23_CHOCOPhlAnSGB_202307` | Bowtie2 index name inside --metaphlan_db |
| `--humann3_nucleotide_db` | `null` | Set in a site profile or with --humann3_nucleotide_db |
| `--humann3_protein_db` | `null` | Set in a site profile or with --humann3_protein_db |
| `--adapters` | `ref/known_adapters.fna` | the adapter FASTA for `--adapter_trim` |
| `--adapter_trim` | `auto` | auto | always | never -- see the adapter_trim note in the workflow below |

## Conda Environment Paths

| Parameter | Default | Description |
|---|---|---|
| `--samtools_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--fastp_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--fastqc_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--minimap2_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--multiqc_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--krakenuniq_bracken_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--consensus_taxa_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--metaphlan4_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--humann3_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--decontam_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--validation_env` | `conda_envs/lock/<env>.linux-64.txt` | the locked conda environment for each tool |
| `--batch_corr_env` | pre-built TSCC env | see *Conda/Mamba* |
| `--antismash_env` | `conda_envs/antismash_env.yml`, `ANTISMASH/` | reserved for the antiSMASH module in `Modules/antismash.nf`, which `main.nf` does not currently call |
| `--batch_corr_method` | `tune` | `tune` (auto-phased Tune_ConQuR) or `vanilla` (default ConQuR settings). Nothing else is accepted: ComBat is not implemented |

## Script Paths

| Parameter | Default | Description |
|---|---|---|
| `--scripts` | `scripts/` in the repo | the helper-script directory |
| `--decontam_script` | `scripts/decontamination.R`, `scripts/2500703_batch_correction_normalization.r` | the two R scripts, overridable for development |
| `--batch_corr_script` | `scripts/decontamination.R`, `scripts/2500703_batch_correction_normalization.r` | the two R scripts, overridable for development |

---
[Back to the README](../README.md)
