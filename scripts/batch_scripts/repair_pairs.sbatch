#!/bin/bash

# Re-pair R1/R2 FASTQ files for nf-core/mag input.
# Our host depletion processes R1 and R2 independently, causing pair mismatch.
# This array job repairs all 520 tumor samples in parallel (50 at a time).

set -e
set -o pipefail

PROJ_DIR="/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy"
OUT_DIR="${PROJ_DIR}/RESULTS/MAPPED_READS_PAIRED"

source ~/.bashrc
conda activate microbiome

cd "${PROJ_DIR}"
mkdir -p "${OUT_DIR}"

# Get the sample for this array task (skip header, 1-indexed)
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" mag_samplesheet.csv | cut -d',' -f1)

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

R1="RESULTS/MAPPED_READS/${SAMPLE}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"
R2="RESULTS/MAPPED_READS/${SAMPLE}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"

echo "=========================================="
echo "Job: ${SLURM_JOB_NAME} (task ${SLURM_ARRAY_TASK_ID})"
echo "Sample: ${SAMPLE}"
echo "Started: $(date)"
echo "=========================================="

if [[ ! -f "${R1}" ]] || [[ ! -f "${R2}" ]]; then
    echo "ERROR: Input files not found for ${SAMPLE}"
    echo "  R1: ${R1} (exists: $(test -f ${R1} && echo yes || echo no))"
    echo "  R2: ${R2} (exists: $(test -f ${R2} && echo yes || echo no))"
    exit 1
fi

python scripts/repair_pairs.py \
    --r1 "${R1}" \
    --r2 "${R2}" \
    --out-dir "${OUT_DIR}"

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
