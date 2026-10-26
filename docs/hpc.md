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

## Launch environment

Nextflow copies the launching shell's environment into every task (and `sbatch` copies it into the
driver). Two things in a login shell can make a task run something other than its pinned tool:

- **PATH order.** Nextflow activates a task's env with `source <base>/bin/activate <env>`, which puts
  the env's `bin` where the conda base's `bin` was, not first. A tool directory earlier on PATH
  (on TSCC, `~/.bashrc`'s `tools/minimap2` 2.24, or `nf-env/bin` with samtools, bgzip, krakenuniq,
  bracken and python) then wins over the env. CRC50 attempt 1 (Job 12380606) ran minimap2 2.24 in
  every mapReads task this way, against a pinned 2.28.
- **Library and module paths.** `LD_LIBRARY_PATH`, `LIBRARY_PATH`, `CPATH`, `PYTHONPATH`,
  `PYTHONHOME`, `PERL5LIB` and `R_LIBS*` from the login shell.

The pipeline guards against both. `process.beforeScript` (`nextflow.config`) unsets those variables
and puts `$(conda info --base)/bin` first before activation, so the env's `bin` comes first after it;
`PYTHONNOUSERSITE=1` keeps `~/.local` packages out. Every process then records the version and path
of each tool it runs, and `--tool_version_check strict` (default) stops the run (exit 78,
`terminate`) if one resolves outside the task's env or differs from its lockfile pin. `warn` records
the problem in `RUN_REPORT.txt` and carries on. The provenance finaliser runs with
`LD_LIBRARY_PATH`, `PYTHONPATH` and `PYTHONHOME` dropped, and its exit status is reported.

A clean launch is still the first line of defence: start the driver from `PATH=/usr/local/bin:/usr/bin:/bin`,
source `<conda base>/etc/profile.d/conda.sh`, add only Java, and call Nextflow by absolute path (see the
`env -i` re-exec at the top of the TSCC production drivers, which keeps only `SLURM_*` and a few basics).

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

**Batch correction** is built by Nextflow from `lock/batch_correction_env.linux-64.txt` like the
other environments. Its two packages that no channel provides, r-cqrreg 1.2.1 and r-conqur
(`ivartb/ConQuR_par` at `--conqur_sha`), come from `conda_envs/local_channel/`; the lockfile names
them under `@LOCAL_CHANNEL@`, which the workflow replaces with that directory's absolute path (so
the checkout path may contain only letters, digits and `. _ / + -`). The
workflow checks the installed ConQuR commit before execution. `--batch_corr_env <prefix>` still
selects a pre-built environment (`scripts/setup_batch_env.sh` builds one from the same lockfile).
See [batch setup and reproducibility](running/batch_setup_and_reproducibility.md).

---
[Back to the README](../README.md)
