# Running CMPipeline on HPC systems

CMPipeline submits every process to the batch scheduler through a Nextflow executor profile. All compute
nodes must be able to read the repository, the inputs, the databases, the conda cache, the work directory
and the output directory.

## General recommendations

- Launch Nextflow from a small **driver job** (2 CPUs, 8 GB, a few days), not a login node or a notebook
  session: when a notebook job ends, a driver running inside it dies mid-run.
- Keep the work directory on a shared high-throughput filesystem, and keep the same work directory and
  parameters when using `-resume`.
- Put site accounts, queues and QOS in a private config (`-c my-site.config`), not in the repository.
- For a production run, execute a **frozen copy** of the code (a read-only checkout or snapshot), so that
  editing the repository while a run is going cannot change what running tasks execute.
- Start with one or two samples before a cohort.

## Example driver job (Slurm)

```bash
#!/usr/bin/env bash
#SBATCH --job-name=cmpipeline-driver
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
# every #SBATCH line directly under the shebang: directives after the first command are ignored

nextflow run /path/to/CMPipeline/main.nf -profile slurm -c my-site.config \
    --sample samples.csv --hg38_db /refs/human-GRC-db.mmi --t2t_phix_db /refs/human-GCA-phix-db.mmi \
    --kraken_db /db/krakenuniq --metaphlan_db /db/metaphlan --run_decontam false -resume
```

On TSCC, `scripts/batch_scripts/main.sbatch` is a ready driver: it takes the sample sheet as its first
argument and passes anything further to Nextflow (see `scripts/batch_scripts/README.md`).

## Resource expectations

Every request scales with the retry attempt, so a task that runs out of memory or time gets twice as much
on its next try. Sized from measured peak memory on TCGA CRC WGS where noted.

| Process | Label | CPUs | Memory | Time |
|---|---|---:|---:|---:|
| BAM/CRAM extraction | `extract_reads` | 12 | 4 GB | 8 h |
| fastp filtering | `filter_reads` | 12 | 8 GB | 10 h |
| FastQC | `fastqc` | 4 | 4 GB | 8 h |
| Host depletion (minimap2) | `mapReads` | 4 | 24 GB | 50 h |
| KrakenUniq + Bracken | `process_high_disk` | 4 | 128 GB | 50 h |
| MetaPhlAn 4, HUMAnN 3, batch correction | `process_high` | 16 | 256 GB | 50 h |
| Consensus, decontamination, merges, MultiQC | `process_medium` | 4 | 16 GB | 30 h |
| Reports and small steps | `process_low` | 2 | 4 GB | 4 h |
| Validation reference build | (fixed) | 16 | 96 GB | 12 h |
| Validation alignment | (scales with the index) | 16 | 1.2 × index + 12 GB | 12 h |

Per-task ceilings come from `--max_cpus`, `--max_memory` and `--max_time`. Under `-profile local` they are
8 CPUs and 32 GB.

## Profiles

- `-profile tscc` uses the TSCC Slurm queue/account and TSCC database defaults.
  Pangenome depletion is **not** among them: that index set is large enough that
  enabling it silently would dominate every run, so pass `--pangenome_db` explicitly
  when you want it.
- `-profile slurm` is a generic Slurm profile; override `process.queue` and
  `process.clusterOptions` in a site overlay.
- `-profile pbs` and `-profile lsf` provide generic PBS/Torque and LSF templates.
- `-profile biowulf` targets NIH Biowulf (the `norm` partition, conservative scheduler
  polling, and `--gres=lscratch:` from `--biowulf_lscratch_gb`).
- `-profile sge` is a generic Sun/Oracle Grid Engine template.
- `-profile local` runs on the current machine. It caps total concurrency at 8 CPUs and
  32 GB deliberately: `process.resourceLimits` bounds each individual task, but only
  `executor.cpus`/`executor.memory` stop a local run swamping a shared login node.

Example TSCC launch:

```bash
nextflow run main.nf -profile tscc --sample samples.csv -resume
```

For another cluster, copy the relevant `conf/*.config` file or provide an overlay:

```bash
nextflow run main.nf -profile slurm -c my-site.config --sample samples.csv
```

## Conda environments

Conda is enabled by default. To use Mamba where available:

```bash
nextflow run main.nf -profile local,mamba --sample samples.csv
```

Built environments are cached in `.conda_cache/` next to the workflow rather than
inside the work directory, so they survive a work-directory cleanup and are shared
between runs. Override with `NXF_CONDA_CACHEDIR` or `conda.cacheDir` in a site profile.

Each tool environment is built from an exact lockfile, `conda_envs/lock/<env>.linux-64.txt`
(`conda list --explicit`), selected by the `--<tool>_env` parameters listed under
*All parameters*. The `.yml` beside each lockfile is its human-readable source, and
`conda_envs/build_envs.sbatch` regenerates the lockfiles. A pre-existing environment such as
`nf-env` may be used to launch Nextflow, but tool environments are still built from the lockfiles.

**Batch correction is the exception: it runs in a pre-built environment**,
`--batch_corr_env` (on TSCC,
`/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/envs/batch_correction_20260923`).
ConQuR is not packaged for conda, so an environment built from the lockfile alone would lack it.
`sbatch conda_envs/build_envs.sbatch` builds that environment and runs
`scripts/install_conqur.R` in it, which installs cqrReg 1.2.1 and the `ivartb/ConQuR_par` fork
of ConQuR at the commit in `--conqur_sha`. Before batch correction starts, the workflow checks
that ConQuR is installed **at that exact commit** and stops with the command to fix it
otherwise. On another site, build the environment the same way and point `--batch_corr_env` at it.

---
[Back to the README](../README.md)
