#!/bin/bash

# ============================================================================
# Run antiSMASH on MEGAHIT assembly contigs or MAG bins
#
# Two modes:
#   contigs (default) - Run on all 520 MEGAHIT assemblies directly
#   mags              - Run on binned MAGs from nf-core/mag
#
# Prerequisites:
#   1. MEGAHIT assemblies exist (RESULTS/MAG/Assembly/MEGAHIT/) or MAG bins
#   2. antiSMASH conda env created:
#        mamba env create -f conda_envs/antismash_env.yml
#        conda activate nf-antismash-env
#        download-antismash-databases
#        antismash --prepare-data
#
# Usage:
#   sbatch scripts/batch_scripts/run_antismash.sbatch                  # contigs, fresh
#   sbatch scripts/batch_scripts/run_antismash.sbatch -resume           # contigs, resume
#   sbatch scripts/batch_scripts/run_antismash.sbatch --mode mags       # MAGs
# ============================================================================

set -e
set -o pipefail

echo "=========================================="
echo "Job: ${SLURM_JOB_NAME}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Started: $(date)"
echo "=========================================="

source ~/.bashrc
conda activate microbiome

PROJ_DIR="/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy"
cd "${PROJ_DIR}"

mkdir -p scripts/batch_scripts/log

# Default to contigs mode; pass through any extra args (--mode mags, -resume, etc.)
nextflow run run_antismash.nf \
    --mode contigs \
    -c conf/antismash.config \
    -work-dir work/antismash \
    "$@" \
    -N l1joseph@ucsd.edu

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
