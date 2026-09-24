#!/bin/bash
#SBATCH --job-name=cmpipeline-resume
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/cmpipeline-resume-%j.out
#SBATCH --error=logs/cmpipeline-resume-%j.err

# Every #SBATCH directive must stay above the first executable line -- Slurm stops
# reading them at that point, and the job silently inherits the 1 h default walltime.
# Validate changes with: sbatch --test-only scripts/batch_scripts/main.sbatch
#
# This job is only the Nextflow head process; the pipeline's own work is submitted
# as separate Slurm jobs by -profile tscc. 4 cpus / 16 GB is for the head process.
#
# Usage:
#   cd <your CMPipeline checkout>
#   mkdir -p logs
#   sbatch scripts/batch_scripts/resume_main.sbatch <sample_sheet.csv> [extra nextflow args]
#
# Identical to main.sbatch -- both pass -resume. Kept as a separate name because the
# README and older notes refer to it.

set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$PWD}"
SAMPLE_SHEET="${1:-$REPO/samples.csv}"
shift || true

# Pinned for reproducibility; the config is verified on 24.10.0 and 26.04.6 (nextflow.config).
export NXF_VER="${NXF_VER:-24.10.0}"
CONDA_BASE="${CONDA_BASE:-/tscc/projects/ps-lalexandrov/amabbasi/anaconda3}"
NEXTFLOW="${NEXTFLOW:-$CONDA_BASE/envs/nf-env/bin/nextflow}"
export JAVA_HOME="${JAVA_HOME:-$CONDA_BASE/envs/nf-env/lib/jvm}"
export PATH="$JAVA_HOME/bin:$PATH"

# No launcher env is activated (audit D09): Nextflow and Java come from the explicit nf-env
# paths above, and every process gets its own conda env from conda_envs/lock/. The old
# `conda activate microbiome` named an env no yml defines and aborted the job under set -e.

cd "$REPO"
echo "Repo:         $REPO"
echo "Sample sheet: $SAMPLE_SHEET"
echo "Nextflow:     $($NEXTFLOW -v 2>&1)"

"$NEXTFLOW" run main.nf \
    -profile tscc \
    --sample "$SAMPLE_SHEET" \
    -resume \
    "$@"
