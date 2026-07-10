#!/bin/bash

# ============================================================================
# Download GTDB-Tk and CheckM2 databases for nf-core/mag
#
# GTDB-Tk database: ~85 GB (compressed ~35 GB download)
# CheckM2 database: ~3 GB
#
# Usage:
#   sbatch scripts/batch_scripts/download_gtdb_checkm2.sbatch
# ============================================================================

set -e
set -o pipefail

echo "=========================================="
echo "Job: ${SLURM_JOB_NAME}"
echo "Started: $(date)"
echo "=========================================="

source ~/.bashrc

DB_BASE="/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/databases"
mkdir -p "${DB_BASE}"

# ── 1. Download GTDB-Tk database ──────────────────────────────────────────
GTDB_DIR="${DB_BASE}/gtdbtk"
if [[ -d "${GTDB_DIR}" ]] && [[ -f "${GTDB_DIR}/metadata/metadata.txt" ]]; then
    echo "GTDB-Tk database already exists at ${GTDB_DIR}, skipping"
else
    echo "Downloading GTDB-Tk database (R220) ..."
    mkdir -p "${GTDB_DIR}"
    cd "${DB_BASE}"

    # GTDB-Tk R220 (latest as of 2024)
    wget -c https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz

    echo "Extracting GTDB-Tk database ..."
    tar -xzf gtdbtk_r220_data.tar.gz -C "${GTDB_DIR}" --strip-components=1
    rm -f gtdbtk_r220_data.tar.gz

    echo "GTDB-Tk database installed at ${GTDB_DIR}"
fi

# ── 2. Download CheckM2 database ──────────────────────────────────────────
CHECKM2_DIR="${DB_BASE}/checkm2"
if [[ -d "${CHECKM2_DIR}" ]] && [[ -f "${CHECKM2_DIR}/uniref100.KO.1.dmnd" ]]; then
    echo "CheckM2 database already exists at ${CHECKM2_DIR}, skipping"
else
    echo "Downloading CheckM2 database ..."
    mkdir -p "${CHECKM2_DIR}"

    # CheckM2 uses a diamond database
    conda activate microbiome 2>/dev/null || true
    pip install checkm2 2>/dev/null || true

    # Download via checkm2 if available, otherwise manual
    if command -v checkm2 &>/dev/null; then
        checkm2 database --download --path "${CHECKM2_DIR}"
    else
        # Manual download
        cd "${CHECKM2_DIR}"
        wget -c https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
        tar -xzf checkm2_database.tar.gz --strip-components=1
        rm -f checkm2_database.tar.gz
    fi

    echo "CheckM2 database installed at ${CHECKM2_DIR}"
fi

echo "=========================================="
echo "Database downloads complete"
echo "GTDB-Tk: ${GTDB_DIR}"
echo "CheckM2: ${CHECKM2_DIR}"
echo "Finished: $(date)"
echo "=========================================="
