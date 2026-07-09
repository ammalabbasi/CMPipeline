#!/bin/bash
#SBATCH -p platinum
#SBATCH -q hcp-ddp302
#SBATCH -A ddp302
#SBATCH --job-name=eac_diff_abund
#SBATCH --output=logs/eac_diff_abund_%j.out
#SBATCH --error=logs/eac_diff_abund_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mail-user=l1joseph@ucsd.edu
#SBATCH --mail-type=END,FAIL

set -e
set -o pipefail

echo "=========================================="
echo "Job: ${SLURM_JOB_NAME}"
echo "Started: $(date)"
echo "=========================================="

source ~/.bashrc

PIPELINE=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy
OTU_GENUS=$PIPELINE/RESULTS/05_BATCH_CORRECTION/shipment_rerun/corrected/ConQuR_tuned.tsv
OTU_SPECIES=$PIPELINE/RESULTS/05_BATCH_CORRECTION/species_shipment/corrected/ConQuR_tuned.tsv
META=$PIPELINE/RESULTS/04_DECONTAMINATION/bracken_run/EAC_GEJ_METADATA.txt
OUTDIR=$PIPELINE/RESULTS/06_DA_ANALYSIS
mkdir -p $OUTDIR/logs

cd $PIPELINE

# ============================================================================
# MaAsLin3 - Genus
# ============================================================================
echo ""
echo "--- MaAsLin3 Genus ($(date)) ---"
conda run -n maaslin-env Rscript scripts/run_maaslin.R \
  --otu $OTU_GENUS \
  --meta $META \
  --prefix $OUTDIR/maaslin/genus_type \
  --formula "Type" \
  --level genus \
  --cores 16
echo "MaAsLin3 genus DONE"

# ============================================================================
# MaAsLin3 - Species
# ============================================================================
echo ""
echo "--- MaAsLin3 Species ($(date)) ---"
conda run -n maaslin-env Rscript scripts/run_maaslin.R \
  --otu $OTU_SPECIES \
  --meta $META \
  --prefix $OUTDIR/maaslin/species_type \
  --formula "Type" \
  --level species \
  --cores 16
echo "MaAsLin3 species DONE"

# ============================================================================
# ANCOMBC2 - Genus
# ============================================================================
echo ""
echo "--- ANCOMBC2 Genus ($(date)) ---"
conda run -n ancombc-env Rscript scripts/run_ancombc.R \
  --otu $OTU_GENUS \
  --meta $META \
  --prefix $OUTDIR/ancombc/genus_type \
  --formula "Type" \
  --level genus
echo "ANCOMBC2 genus DONE"

# ============================================================================
# ANCOMBC2 - Species
# ============================================================================
echo ""
echo "--- ANCOMBC2 Species ($(date)) ---"
conda run -n ancombc-env Rscript scripts/run_ancombc.R \
  --otu $OTU_SPECIES \
  --meta $META \
  --prefix $OUTDIR/ancombc/species_type \
  --formula "Type" \
  --level species
echo "ANCOMBC2 species DONE"

# ============================================================================
# BIRDMAn NB - Genus  (slowest step — ~4-8h for full dataset)
# ============================================================================
echo ""
echo "--- BIRDMAn NB Genus ($(date)) ---"
conda run -n q2-birdman python scripts/run_birdman.py \
  --otu $OTU_GENUS \
  --meta $META \
  --prefix $OUTDIR/birdman/genus_type \
  --formula "C(Type, Treatment('Normal'))" \
  --level genus \
  --num_iter 500
echo "BIRDMAn genus DONE"

# ============================================================================
# BIRDMAn NB - Species
# ============================================================================
echo ""
echo "--- BIRDMAn NB Species ($(date)) ---"
conda run -n q2-birdman python scripts/run_birdman.py \
  --otu $OTU_SPECIES \
  --meta $META \
  --prefix $OUTDIR/birdman/species_type \
  --formula "C(Type, Treatment('Normal'))" \
  --level species \
  --num_iter 500
echo "BIRDMAn species DONE"

# ============================================================================
# Concordance - Genus
# ============================================================================
echo ""
echo "--- DA Concordance Genus ($(date)) ---"
conda run -n da_concordance python scripts/run_da_concordance.py \
  --maaslin  $OUTDIR/maaslin/genus_type_maaslin_full.tsv \
  --ancombc  $OUTDIR/ancombc/genus_type_ancombc2_full.tsv \
  --birdman  $OUTDIR/birdman/genus_type_birdman_results.tsv \
  --level genus \
  --output   $OUTDIR/concordance/da_concordance_genus.tsv
echo "Concordance genus DONE"

# ============================================================================
# Concordance - Species
# ============================================================================
echo ""
echo "--- DA Concordance Species ($(date)) ---"
conda run -n da_concordance python scripts/run_da_concordance.py \
  --maaslin  $OUTDIR/maaslin/species_type_maaslin_full.tsv \
  --ancombc  $OUTDIR/ancombc/species_type_ancombc2_full.tsv \
  --birdman  $OUTDIR/birdman/species_type_birdman_results.tsv \
  --level species \
  --output   $OUTDIR/concordance/da_concordance_species.tsv
echo "Concordance species DONE"

echo ""
echo "=========================================="
echo "All DA analysis complete: $(date)"
echo "Results: $OUTDIR"
echo "=========================================="
