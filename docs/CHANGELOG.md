# Changes in this version

Not yet released; these are the changes since the last tagged state of `main`.

**Correctness fixes that change results** — read these before comparing against older output.

- **June-Young's known-contaminant list is back, as a flag (`--known_contaminants`, `ref/known_contaminants.tsv`).**
  His decontam script carried 95 species (`blacklist_species`; mostly reagent genera such as *Acinetobacter*,
  *Methylobacterium*, *Sphingomonas*, *Ralstonia*) and used them to label plots. The 24 Sep rewrite of
  `decontamination.R` dropped it. It is restored verbatim and, at species rank, flags listed taxa
  (`known_contaminant`) and reports each sample's share of reads on the list (`pct_known_contaminant`,
  Ludmil's reagent check). It removes nothing; decontam's calls are unchanged.
- **Human and PhiX are removed before decontam again (`--decontam_remove_taxa`).** June-Young's script removed
  *Homo sapiens* and PhiX first (`hardcoded_removal_list`); the 24 Sep rewrite dropped that. Default
  `Homo sapiens;Sinsheimervirus phiX174` (KrakenUniq's names; genera `Homo` / `Sinsheimervirus` at genus rank),
  exact-name match instead of his substring match. **Decontam input, and so its calls, can change.**

- **Validation reference is ANI-curated; new `validated_group` status (decision 36).** Genomes RefSeq filed
  under the wrong species move to the species they belong to (4 *Fusobacterium* "nucleatum" genomes are
  *animalis*/*vincentii*; 4 *Bacteroides* and 4 *Escherichia* `sp.` genomes are named species). Species
  whose anchors meet at 92% ANI form a fallback group (*F. nucleatum* group: *nucleatum*, *animalis*,
  *polymorphum*, *vincentii*, *hwasookii*, *canifelinum*). A cell that is `insufficient` only because its
  reads are shared, and whose group validates, is `validated_group` (count kept). The unnamed 18-genome
  cluster inside *F. nucleatum* stays there. **A new reference is built (cache key v5)**, and validation
  statuses change, mostly for *Fusobacterium*. New stats columns `ani_group`, `group_status`,
  `group_fail_reason`; `validation_mask_report.tsv` gains `n_validated_group`.

- **Reads extracted from BAM/CRAM are adapter-trimmed again (`--adapter_trim auto`).** Since 16 Sep 2026
  `auto` passed no adapter FASTA for them, assuming aligned reads were trimmed. GDC did not trim TCGA, and
  bwa soft-clips adapter read-through, so it reached KrakenUniq as *synthetic construct*. `auto` now uses a
  5-sequence Illumina list (`--adapters_core`, `ref/illumina_core_adapters.fna`) for them; FASTQ keeps the
  234-sequence list. Synthetic test: all planted adapters found at 17 s vs 645 s for the full list.
  **Read counts after fastp, and synthetic-construct counts, change for BAM/CRAM input.**

- **QC summary: human and synthetic construct reported apart from microbial reads.** New
  `kraken_classified_reads` (the root count the MetaPhlAn gate uses, unchanged), `kraken_human_reads`,
  `kraken_synthetic_reads`, and Bracken species-level `bracken_human_reads`, `bracken_synthetic_reads`,
  `pct_bracken_human`, `pct_bracken_synthetic`. **`bracken_microbial_reads` now excludes human and
  synthetic construct**; before, it was the root count. `low_as_share_of_bracken` now divides by
  `kraken_classified_reads` (same value as before). New `adapter_list`; `adapter_scan_applied` is read from
  fastp's command line, because fastp writes `adapter_cutting` whenever it trims, auto-detection included (0.24 too). The old
  column was therefore wrong for BAM/CRAM input: `true` meant "auto-detection trimmed something".

- **Tool update (2026-10-26):** minimap2 2.28 → 2.31, samtools → 1.24, fastp 0.24.0 → 1.3.7, Bracken 2.9 →
  3.1p1, bowtie2 (MetaPhlAn env) → 2.5.5; the Kraken/Bracken env is reduced to the tools the task runs.
  These are the versions pksProfiler validated. **fastp 1.x can change read counts** (poly-G trimming and
  adapter matching). The R/ConQuR/decontam stack, MetaPhlAn 4.1.1 + its database, HUMAnN, skani and the
  KrakenUniq database stay frozen, with reasons in `docs/tool_versions.md`.

- **Bracken now uses each sample's read length (`--bracken_read_length auto`, new default).** Bracken
  re-estimates species abundance with a k-mer distribution built for one read length. The default was a
  fixed 50, applied to ~150 bp TCGA WGS reads although the database ships 50–300 bp distributions. `auto`
  takes the median read length of the classified reads and the closest available distribution, and writes
  the choice to `<sample>.bracken_read_length.tsv`. **Species-level Bracken results change**; an integer
  restores a fixed length.

- **Bracken failed for every sample since `--bracken_read_length auto` (`1cf242b`).** The
  `kmer_distrib` sed in `Modules/Bracken.nf` had single backslashes before `(`, `)`, `.` and `1`
  inside the Groovy script string. Nextflow's v1 parser (24.10/25.x default) refused to compile
  `main.nf`; the v2 parser (26.04) compiled it but shifted every later interpolation of the Bracken
  script onto the wrong value, so `kreport2mpa.py` was called from the sample name and every task
  failed. The backslashes are doubled; `tests/test_bracken_read_length.py` now scans every module
  script for invalid Groovy escapes. No run produced results with the broken line.

- **Alignment validation is also run with `--frag=no`.** `alignForValidation` streamed interleaved mates
  into `minimap2 -ax sr` in fragment mode, the same input-route dependence fixed for host depletion.
  Validation calls may change; re-score.

- **Batch correction no longer stops a whole rank when one prevalence bin sits in a single batch
  (audit K04).** `Tune_ConQuR` scored each bin's uncorrected table with `adonis(tab ~ batch)` without
  the single-batch guard it uses for corrected tables, so one batch-confined rare taxon failed the rank
  with "contrasts can be applied only to factors with 2 or more levels" (CRC50 Test 3, genus). The
  pipeline now calls a guarded copy (`guarded_tune_conqur`, `scripts/batch_contract.R`). Such a bin is
  scored on its corrected candidates only. If none can be scored, its taxa stay uncorrected and
  `batch_correction_methods.tsv` records `original_selected`. Runs that completed before are
  unchanged: the six synthetic fixtures that failed now complete, and a passing fixture gives
  byte-identical output. The packaged ConQuR is not modified; the patch stops if the pinned code differs.

- **Host depletion no longer depends on input route or record order (`--host_depletion_frag no`, default).**
  minimap2's `-ax sr` preset pairs adjacent mate-named reads (fragment mode). Reads extracted from a BAM/CRAM
  have mates adjacent, paired FASTQ has all R1 then all R2, so the same reads were depleted differently by
  route. Every host-depletion call now runs with `--frag=no`. Tested on the real GRCh38/T2T indexes: only
  human-like reads changed (BAM-order runs had removed more of them); microbial reads were identical. BAM/CRAM
  cohorts processed before this change should be re-run from host depletion. `--host_depletion_frag auto`
  restores the old behaviour.

- **BAM/CRAM extraction missed low-score mapped reads in `bwa mem -T 0` alignments (Ludmil 1 Oct,
  L8; ported from pksProfiler E1).** GDC aligned TCGA with `bwa mem -T 0`, so microbial reads with a
  chance short human hit were reported mapped (AS < 30) and FLAG 4 extraction never saw them.
  `--extract_low_score_mapped auto|on|off` (new, default `auto`) also extracts primary mapped records
  with `AS:i` < 30; `auto` is on only when the input's `@PG` runs `bwa mem`/`bwa-mem2 mem` with
  `-T` < 30, and off for STAR, minimap2, `bwa aln`, default-`-T` bwa mem and inputs without `@PG`.
  Secondary/supplementary records and records without `AS` are never added. **Results change for
  TCGA GDC BAMs** (more reads reach classification); other inputs, and `off`, extract exactly as
  before. Ten new QC summary columns (48 in all) record the upstream aligner, version and `-T`, the
  mode and why, and the rescued reads extracted / after host depletion / classified; a new
  `lowAsSurvival` task counts the last two by read name. `RUN_REPORT.txt` gains a counts-only
  *LOW-SCORE RESCUE* section. Per-taxon attribution through Bracken is not measured (see
  docs/running/classification.md).

- **Alignment validation judged spread on reads shared with relatives (audit 2026-09-30, V01/V02).**
  The bins test counted every supporting read, including reads tied with another species. For a
  species present at low abundance next to an abundant relative, those shared reads pile into the
  regions the two genomes have in common, so the species failed `bins_ratio` although it was there
  (22 of 42 such cases in a synthetic benchmark on public genomes), while absent relatives could
  validate (18 of 58). Bins and the reads-on-genome gate are now computed on reads unique to the
  species (`--validation_bins_reads unique`, new, default), and `--validation_min_unique_frac` is 0.05
  (was 0.02). Benchmark: minors 29 validated / 4 failed, decoys 1 of 58 validated, uniform species
  129 of 131. `--validation_bins_reads all --validation_min_unique_frac 0.02` restores the old rule.
  New stats columns `bins_reads`, `reads_for_bins`. Identity 0.98 is unchanged: on CRC50 the identity
  failures sat at 4–10 mismatches per 151 bp, validated cells at 0–3.

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

**phiX: the "T2T + PhiX" labels were wrong; the index holds no phiX** (port review #19). The second
host-depletion index, `human-GCA-phix-db.mmi` (`--t2t_phix_db`), is T2T-CHM13v2.0 alone: 25 sequences,
the published CHM13v2.0 total to the base (read from the index by pksProfiler, 2026-10-22). In the code,
`scripts/download_references.sh` concatenates `ref/NC_001422.fna` into the index FASTA but never downloads
it, and has no `set -e`. The cohort report said "Removed: T2T + PhiX" and the README "T2T-CHM13+PhiX".
- Labels now say T2T-CHM13; `docs/host_depletion_references.md` documents both indexes, the mechanism, and
  what it means: **phiX spike-in reads are not removed and count in KrakenUniq/Bracken totals and library
  sizes**; viruses are dropped at consensus, so consensus and downstream tables carry no phiX taxon (with
  `--skip_consensus`, the Bracken tables that go on to decontamination do).
- Unchanged on purpose: the index, `--t2t_phix_db`, and the columns `reads_after_t2t_phix` /
  `pct_removed_t2t_phix` (rebuilding mid-cohort would change depletion for samples already run). No
  result changes.

**QC floor: a sample that cannot be assessed is `not_assessable`, never a microbial negative**
(pksProfiler COHORT_CHECKS A3, a proposal there; port review #34). A BAM/CRAM pre-filtered to mapped
reads gives an empty profile that decontam and ConQuR would read as "no microbes" and that would sit in
the low-biomass tail.
- `extractReads` records `input_qc_floor` in `<sample>.extract_counts.tsv`: `unmapped_filtered_upstream`
  when the input's `@PG` chain shows its unmapped reads were removed upstream (`input_unmapped_filter`:
  `samtools view -F 4`, `-f 2`, `-q` > 0, a region; a sambamba filter), `no_unmapped_records` when the
  BAM/CRAM has no FLAG 4 record (the low-score rescue does not lift this), else `pass`.
- `--min_reads_after_depletion` (new, default `0` = off): fewer host-depleted reads than this is
  `below_min_reads_after_depletion`. No calibrated value exists yet, so it is off.
- `--qc_floor_action exclude|flag` (new, default `exclude`): `exclude` stops a not-assessable sample at
  extraction (or before classification, for the depth floor), so it is in no cohort table; `flag` profiles
  it and only marks it.
- QC summary: `terminal_state` `not_assessable_<reason>`, `sample_status` `not_assessable`, new column
  `qc_floor` (column count then: 64). `RUN_REPORT.txt` *SAMPLE OUTCOMES* states the count and that they are not
  negatives; the cohort report adds a *Not assessable* tile.
- **Changes which samples are reported, not any count**: a BAM/CRAM with zero unmapped records, or with
  an upstream unmapped filter in its `@PG`, no longer reaches the cohort tables by default (zero-record
  inputs were already dropped at fastp; filtered inputs with rescued low-AS reads were profiled before).
- The channel helpers live in `lib/QcFloor.groovy`: with them inline, `main.nf`'s entry workflow
  exceeded the 65,535-character string Nextflow 24.10's parser allows. `tests/test_repo_hygiene.py` now
  checks that limit (review #38); the block is at ~64,700 characters, so new workflow logic belongs in `lib/`.

**Input provenance for BAM/CRAM: reference fingerprint and aligner chain** (audit X13, Ludmil E3;
pksProfiler `70cb629`). Nothing recorded which human build each BAM/CRAM was aligned to, or how, so a
results tree mixing TCGA (GDC GRCh38.d1.vd1: hs38d1 decoys, chrEBV, ~200 viral contigs), Mutographs and
Hartwig inputs could not say that their unmapped reads were defined against different contig sets.
- `extractReads` writes `<sample>.input_provenance.tsv` from the header it already reads for the
  low-score rescue decision (`low_score_rescue.py decide --provenance`, `scripts/header_provenance.py`):
  the header is read once per sample.
- Twelve QC summary columns (63 in all): `input_ref_fingerprint` (md5 of ordered `@SQ` SN:LN),
  `input_ref_m5_fingerprint`, `input_ref_contigs`, `input_ref_build`, ALT / decoy / HLA / viral contig
  counts, `input_ref_has_ebv`, `input_pg_chain` (program:version), `input_bwa_options` (letters and
  numbers only, e.g. `-T 0`) and `input_unmapped_filter` (an upstream `samtools view -F 4`-type step).
- `RUN_REPORT.txt` gains a counts-only *INPUT ALIGNMENTS* section (`outputs.input_alignments` in the
  JSON): e.g. "50 input(s): 1 reference fingerprint(s); 1 aligner configuration(s)" and
  "50  bwa 0.7.15-r1140 -T 0".
- Never recorded: the `@PG` CL itself, `@RG` text, read-group/sample barcodes, paths. No result changes;
  `extractReads` re-runs on `-resume` (its script and helper digests changed).

**Database manifest: every run records which database releases it read** (audit X14; pksProfiler
`c08098d`). A path is not a version: the KrakenUniq, MetaPhlAn, HUMAnN and host-index directories could
be rebuilt in place and every run would print the same line.
- `conf/databases.tsv` freezes `hg38_db`, `t2t_phix_db`, `pangenome_db`, `kraken_db`,
  `bracken_kmer_distrib`, `metaphlan_db` and both HUMAnN databases: release label, path, size, file
  count, newest mtime, a `sha256:listing` fingerprint (relative path + size of every file; identical to
  pksProfiler's method and values for the shared KrakenUniq and host-index files) and a key check
  (`sha256:files` of small release-defining files, or `mmi:dictionary` for minimap2 indexes).
- `scripts/database_versions.py` recomputes them at the end of each run for the databases the requested
  stages used; `RUN_REPORT.txt` gains *DATABASE VERSIONS* (release, fingerprint, `frozen` /
  `CHANGED: listing` / `CHANGED: key` / `not_in_manifest` / `unused`) and the JSON record `databases`.
  Reads metadata, a few small key files (largest: the 98 MB MetaPhlAn `.pkl`, 125 MB KrakenUniq
  `taxDB`) and the first kB of each `.mmi`; a few seconds per run. A record, not a gate. No result changes.
- The `t2t_phix_db` row states that the index contains no phiX (see the phiX entry above).

**Per-sample failure policy: every sample ends in one reported state** (audit X05/X06; pksProfiler
`dc5d703`, `e126439`).
- Alignment validation batches only samples that are columns of the species table. Under
  `--sample_failure_strategy ignore`, a sample whose Bracken task failed was batched anyway, the scorer
  exited 1 on it, and the strict validation step ended the run; such a sample is now logged and skipped.
- The QC summary gains `metaphlan_status`, `humann_status` and `terminal_state` (51 in all at that point):
  `completed_profiled` / `completed_no_call` / `completed`, `failed_at_<stage>` or `excluded_<reason>`.
  MetaPhlAn and HUMAnN failures, silent before, are named per sample. `main.nf` passes the requested
  stages (`--stages`), each sample's input type, and which HUMAnN profiles finished.
- `RUN_REPORT.txt` gains a counts-only *SAMPLE OUTCOMES* section (`outputs.sample_outcomes` in the JSON).
- `tests/test_sample_failure_policy.py` discovers every per-sample process from the source (first input
  `tuple val(sampleID)`) and fails if one lacks label `per_sample`, or a cohort step carries it.
- **The default stays `finish`** (conservative; the audit asked for X05 to be fixed before deciding).
  `ignore` is now safe and recommended for cohort-scale runs. No result changes: per-sample outputs and
  cohort tables are unchanged for runs where every sample succeeds; validation batch composition is
  unchanged unless a sample is missing from the species table.

**Tool environment guard and runtime tool versions** (audit K06/K07/X04; ported from pksProfiler
F16, `26e196d`, `2bc7ca9`, `2c0f8c8`). CRC50 attempt 1 (Job 12380606) ran `~/.bashrc`'s
`tools/minimap2` 2.24 in all 50 `mapReads` tasks against a pinned 2.28, and nothing in the run record
said so: `RUN_REPORT.txt` listed what the lockfiles specify.
- `process.beforeScript` puts the conda base `bin` first before activation, so each task's env `bin`
  is first after it, and unsets `LD_LIBRARY_PATH`, `LIBRARY_PATH`, `CPATH`, `PYTHONHOME`, `PERL5LIB`
  and `R_LIBS*` leaked from the launching shell (`PYTHONNOUSERSITE`/`PYTHONPATH` were already set).
- Every process that runs a tool records each tool's version and resolved path from inside the task,
  after activation (`scripts/record_tool_versions.sh`, `lib/ToolVersions.groovy`), into
  `<work dir>/.tool_versions.tsv`, and checks it against the task's own lockfile.
  **`--tool_version_check strict` (new, default)** stops the run (exit 78 -> `terminate`, also under
  `--sample_failure_strategy ignore`) when a tool resolves outside the env, is missing, or is off its
  pin; `warn` records and carries on. A version read is retried (stdout + stderr, timeout, backoff).
- `RUN_REPORT.txt` gains *TOOL VERSIONS (as run)* and *TOOL VERSION MISMATCH*; the JSON record gains
  `environment.tool_versions`. The finaliser also flags tasks that left no record (`MISSING`) and
  binaries that disagree with their env's conda-meta (`INSTALLED_MISMATCH`), and exits 1 on any
  problem in strict mode.
- The finaliser now runs with `LD_LIBRARY_PATH`/`PYTHONPATH`/`PYTHONHOME` dropped (the login
  `LD_LIBRARY_PATH` broke its `import hashlib`), and its exit status is reported instead of ignored.
- `trace.txt` adds `attempt`, `cpus`, `memory` and `workdir`.
- MetaPhlAn's `[DEBUG]` PATH/`which`/`--version` echo is replaced by the record.
- Results are unchanged for a run that already used the pinned tools. **Every task's script changed
  (one added line), so `-resume` across this change re-runs every task**; switching
  `--tool_version_check` also changes the scripts and so re-runs them.

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
- MultiQC report, a per-sample QC summary (then 38 fields), and a self-contained cohort report.
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
- Required databases are checked per stage: GRCh38 and T2T-CHM13 (`--t2t_phix_db`) for host depletion,
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
- The batch-correction environment is built by Nextflow from its lockfile like the others.
  ConQuR (the pinned `ivartb/ConQuR_par` fork) and cqrReg 1.2.1, which no conda channel provides,
  are built from `conda_envs/recipes/` and shipped in `conda_envs/local_channel/`; no pre-built
  prefix is needed. The workflow checks the installed commit against `--conqur_sha` before it starts.
  Their upstream sources and license texts are committed alongside
  (`conda_envs/local_channel/NOTICE.md`). The checkout path may contain only letters, digits and
  `. _ / + -`, since it is written into the lockfile as a `file://` URL.

**Infrastructure**
- Resource ceilings use `process.resourceLimits` rather than a helper function in
  `nextflow.config`, so every configuration file parses under Nextflow 24.10 and 25+.
- `scripts/batch_scripts/main.sbatch` and `resume_main.sbatch` carry Slurm directives
  again, run in the directory they are submitted from, and take the sample sheet as an
  argument. The remaining scripts in that directory are documented in
  `scripts/batch_scripts/README.md`; they are not submittable as they stand.

---
[Back to the README](../README.md)
