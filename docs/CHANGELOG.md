# Changes in this version

Not yet released; these are the changes since the last tagged state of `main`.

**Correctness fixes that change results** — read these before comparing against older output.

- **Read extraction was dropping reads and writing sequences into task logs.**
  `extract_reads.nf` called `samtools fastq -N -o FILE`. With `-o`, samtools writes only
  READ1/READ2 to the file and sends **category-0 records** — reads flagged neither READ1 nor
  READ2 — to stdout, which in a Nextflow task is `.command.out`. Those reads were missing
  from the extracted FASTQ, and read sequences were written into every task log. Every
  BAM/CRAM input was affected; pre-extracted FASTQ was not. All categories are now streamed
  into `bgzip`. **Extracted counts can only go up, so any earlier result is a lower bound,
  and existing `work/` directories may contain read sequences.**
- **Duplicate sample IDs passed straight through ingestion.** Two rows with the same
  `patient` value — or IDs differing only by surrounding whitespace, which is stripped —
  collided in `publishDir` and in every merged table, one sample silently overwriting the
  other. Both are now rejected by name, as is a sample sheet with no data rows.
- **`--decontam_min_prevalence` is now 0.05**, matching the published method. It was 0.02.
- **Decontamination is a new implementation** (`scripts/decontamination.R`): blood-derived
  vs solid-tissue normals are distinguished through `--control_role`, contaminants found in
  ≥ `--decontam_min_batches` batches are removed globally rather than zeroed per batch, and
  a taxon that is present but statistically unevaluable is recorded as `unevaluable`, not
  silently as clean. By default (`--decontam_unevaluable_policy skip`) the run continues and the
  across-batch rule decides from the batches that could test it; `error` stops the run instead.
  All-zero samples are kept as evidence of absence (`--decontam_zero_total_policy keep`) unless
  their library was empty.
- **Decontamination from a consensus table could never have run.**
  `consensus_taxa(...).set{}` captured all three process outputs as a multi-channel object,
  so every downstream operator failed. Fixed, and the table is now taken from the channel
  rather than re-read from `publishDir`, which was also a race.

**Boolean flags** — on Nextflow 26, `--flag false` arrives as the String `"false"`, which is
truthy in Groovy, so every boolean passed as `false` would silently have been **on**,
including `--run_decontam false`. Each flag is now resolved once through a single coercion
and rejected outright if it is not a boolean.

**Reproducibility**
- Runs on Nextflow 24.10.0 **and** 26.04.6; no `NXF_VER` pin. `manifest.nextflowVersion` is
  `>=24.10.0`.
- Every run writes a provenance record and `RUN_REPORT.txt`, stating the commit and whether
  the working tree was dirty.
- `trace`, `report`, `timeline` and `dag` are on by default, under the launch directory.
- Conda environment paths are `${projectDir}`-relative; previously they resolved against the
  launch directory, so any run started outside the repository failed to find every
  environment.

**Reporting**
- MultiQC report, a 36-column per-sample QC summary, and a self-contained cohort report.
  `mapReads` now records read retention at each depletion stage, and `filterReads` keeps the
  fastp JSON it previously produced and discarded.

**Input handling**
- BAM, CRAM and FASTQ input through one sample sheet, with `--input_data_type auto`
  detecting each row independently, so one sheet may mix types.
- CRAM sequence names, lengths and MD5s are checked against the reference FASTA before
  decoding, so a mismatched reference is reported instead of producing garbage reads.
- Single-end and paired-end FASTQ; `fastq_1`/`fastq_2` or `r1`/`r2`, second mate optional.
- `--pangenome_db` is optional, and accepts one `.mmi` or a directory of them.

**Failing early instead of late**
- Every path parameter is resolved and checked before the first task is submitted. A flag
  passed without a value (`--pangenome_db` on its own, which Groovy turns into the Boolean
  `true`) is named as such rather than reaching minimap2 hours later.
- Required databases are checked per stage: GRCh38 and T2T+PhiX for host depletion,
  KrakenUniq and MetaPhlAn for classification, both HUMAnN3 databases when it is enabled.

**Speed and cost**
- `--adapter_trim auto|always|never`. `auto` skips the 234-sequence adapter FASTA scan for
  reads extracted from an alignment, which are already trimmed: 2 h 50 m versus 3 m 15 s
  on one aligned CRAM, for zero reads trimmed.
- `-profile tscc` no longer enables pangenome depletion by default.
- Conda environments are cached outside the work directory, so a run no longer rebuilds
  them because it was given a different name.

**Robustness**
- A sample with no reads after filtering is logged and dropped; the run continues. fastp
  exits 0 on empty input while leaving an invalid 0-byte file, which previously failed the
  next step and took the whole run with it.
- HUMAnN3 receives each sample's MetaPhlAn profile as an input instead of reading the
  published output directory while MetaPhlAn was still writing to it.
- `--metaphlan_index` replaces the index name that was hardcoded in two modules.
- Batch correction runs in a pre-built environment (`--batch_corr_env`) built by
  `conda_envs/build_envs.sbatch`, which installs the pinned ConQuR fork with
  `scripts/install_conqur.R`; the conda YAML cannot, as ConQuR is not packaged for conda. The
  workflow checks the installed commit against `--conqur_sha` before it starts.

**Infrastructure**
- Resource ceilings use `process.resourceLimits` rather than a helper function in
  `nextflow.config`, so every configuration file parses under Nextflow 24.10 and 25+.
- `scripts/batch_scripts/main.sbatch` and `resume_main.sbatch` carry Slurm directives
  again, run in the directory they are submitted from, and take the sample sheet as an
  argument. The remaining scripts in that directory are documented in
  `scripts/batch_scripts/README.md`; they are not submittable as they stand.

---
[Back to the README](../README.md)
