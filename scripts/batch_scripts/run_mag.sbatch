#!/bin/bash

# ============================================================================
# Run nf-core/mag on tumor samples (host-depleted reads)
#
# This is the Nextflow head process — it submits individual SLURM jobs
# for each pipeline step via the conf/mag.config executor settings.
#
# Prerequisites:
#   1. Generate samplesheet: python scripts/create_mag_samplesheet.py
#   2. Download databases (one-time):
#      - GTDB-Tk: see scripts/batch_scripts/download_gtdb.sbatch
#      - CheckM2: see scripts/batch_scripts/download_checkm2.sbatch
#   3. Singularity module loaded (auto-loaded below)
#
# Usage:
#   sbatch scripts/batch_scripts/run_mag.sbatch           # fresh run
#   sbatch scripts/batch_scripts/run_mag.sbatch --resume   # resume
# ============================================================================

set -e
set -o pipefail

echo "=========================================="
echo "Job: ${SLURM_JOB_NAME}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Started: $(date)"
echo "=========================================="

# Activate conda environment with Nextflow
source ~/.bashrc
conda activate microbiome

# Load singularity
module load singularitypro/3.11

# Navigate to project directory
PROJ_DIR="/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy"
cd "${PROJ_DIR}"

# Create directories
mkdir -p scripts/batch_scripts/log
mkdir -p /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/singularity_cache

# Database paths
GTDB_DB="/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/databases/gtdbtk"
CHECKM2_DB="/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/databases/checkm2/uniref100.KO.1.dmnd"

# Use paired samplesheet (reads repaired after independent host depletion)
SAMPLESHEET="mag_samplesheet_paired.csv"

# Generate paired samplesheet if it doesn't exist
if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "Generating paired samplesheet..."
    python scripts/update_mag_samplesheet_paired.py
fi

if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "ERROR: ${SAMPLESHEET} not found. Run repair_pairs.sbatch first."
    exit 1
fi

SAMPLE_COUNT=$(tail -n +2 "${SAMPLESHEET}" | wc -l)
echo "Processing ${SAMPLE_COUNT} tumor samples (paired reads)"

# Handle resume flag
RESUME_FLAG=""
if [[ "$1" == "--resume" ]] || [[ "$1" == "-resume" ]]; then
    RESUME_FLAG="-resume"
    echo "Resuming previous run"
fi

# Run nf-core/mag (matching Michael Iter's configuration)
nextflow run nf-core/mag \
    -profile singularity \
    -c conf/mag.config \
    --input "${SAMPLESHEET}" \
    --outdir RESULTS/MAG \
    --skip_shortread_qc \
    --skip_clipping \
    --skip_prodigal \
    --skip_prokka \
    --skip_spades \
    --min_contig_size 2000 \
    --gtdb_db "${GTDB_DB}" \
    --run_checkm2 \
    --checkm2_db "${CHECKM2_DB}" \
    --busco_db_lineage auto_prok \
    -qs 20 \
    -work-dir work/mag \
    ${RESUME_FLAG} \
    -N l1joseph@ucsd.edu

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
