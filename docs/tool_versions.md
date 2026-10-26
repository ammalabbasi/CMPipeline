# Tool versions: what is pinned, what was updated, what is deliberately frozen

Every conda environment runs from an exact lockfile (`conda_envs/lock/<env>.linux-64.txt`,
`conda list --explicit --md5`). The `.yml` files are the human-readable source the lockfiles are solved from.

## Updated 2026-10-26 (branch `feat/tool-bump`)

| Env | Tool | Before | Now | Why / what can change |
|---|---|---|---|---|
| minimap2_env | minimap2 | 2.28 | 2.31 | Latest. Same flags. Host depletion. |
| minimap2_env, samtools_env, metaphlan4_env | samtools | 1.21 / 1.22 | 1.24 | Latest. Same subcommands. |
| validation_env | minimap2 | 2.28 | 2.31 | Same as host depletion, so both alignments use one version. |
| fastp_env | fastp | 0.24.0 | 1.3.7 | Version pksProfiler validated (1.4.0 exists, not yet validated). **Can move read counts:** fastp 1.x enables poly-G trimming for more read-name prefixes and allows one indel in adapter matching. The JSON layout read by `filter_reads.nf` is unchanged. |
| krakenUniq_bracken_env | bracken | 2.9 | 3.1p1 | Same flags (`-d -i -o -w -r -l -t`). Prints `Bracken v3.0.1` (upstream did not update VERSION). In pksProfiler, the abundance and mpa tables were identical; the per-node column of the `-w` report changed. The env was also reduced from a 134-package export to what the task runs (python 3.12, krakenuniq, bracken, gawk, gzip, coreutils). |
| metaphlan4_env | bowtie2 | 2.5.4 | 2.5.5 | Patch release; reads the same mpa_vJun23 index. |

These match the versions pksProfiler adopted on 12 Oct 2026, after a 10-sample TCGA regression with 0 tier changes.
**CMPipeline needs its own before/after regression on CRC50 before tagging** (see `audits/2026-10-26_toolbump/`).

## Deliberately frozen

| Env | Tool(s) | Pinned | Why frozen |
|---|---|---|---|
| batch_correction_env | R 4.3.3, vegan 2.6-4, ConQuR (local build, ivartb/ConQuR_par@ff233085), cqrReg 1.2.1 | as is | ConQuR calls `vegan::adonis`, which later vegan releases removed. ConQuR and cqrReg are built locally for R 4.3 (`conda_envs/local_channel`, `params.conqur_sha`). Changing R or vegan breaks batch correction. |
| decontam_env | R 4.3, decontam 1.22.0, phyloseq 1.46.0 | as is | Bioconductor versions are tied to the R version. Kept on the same R as batch correction, and validated by the decontam tests. Newer decontam (1.30) needs R ≥ 4.5; revisit only together with the R stack. |
| metaphlan4_env, humann3_env | MetaPhlAn 4.1.1 + mpa_vJun23 database | as is | A newer MetaPhlAn (4.2.x) normally comes with a new marker database. That changes consensus calls and needs its own decision, together with the KrakenUniq database version. |
| humann3_env | HUMAnN 3.9, bowtie2 2.5.4, DIAMOND 2.1.12 | as is | HUMAnN 3.9 is the latest release. DIAMOND is left at the version that built/reads the HUMAnN protein database, because `.dmnd` compatibility across DIAMOND releases has not been tested here. HUMAnN is optional in CMPipeline. |
| validation_env | skani 0.3.1 | as is | skani builds the 99.5% dedup and the ANI clusters of the validation reference. A new version can change those clusters and so the reference itself; update only with a reference rebuild. |
| krakenUniq_bracken_env | KrakenUniq 1.0.4 | latest | Already the latest release. |
| — | KrakenUniq database 2023-08-08, MetaPhlAn mpa_vJun23 | as is | Database updates change results more than tool updates. Separate decision. |

## Rules
- Change a pin only together with its lockfile (`conda_envs/build_envs.sbatch`, or the dated build script in the audit folder).
- A tool update that can move counts needs a before/after regression on real samples before a release is tagged.
- When the run-time tool-version guard is enabled, a tool that prints a different version string than its conda
  version (Bracken 3.1p1 → `v3.0.1`) must be listed as a known alias, or the strict check stops the run.
