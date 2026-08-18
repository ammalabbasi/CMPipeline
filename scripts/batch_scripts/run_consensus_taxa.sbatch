#!/bin/bash

set -e
set -o pipefail

echo "=========================================="
echo "Job: ${SLURM_JOB_NAME}"
echo "Consensus Taxa (MetaPhlAn x Bracken)"
echo "Started: $(date)"
echo "=========================================="

source ~/.bashrc
conda activate microbiome

BASE=/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy

python ${BASE}/scripts/compute_consensus_taxa.py \
    --metaphlan ${BASE}/RESULTS/METAPHLAN4/merged_abundance_table_genus.txt \
    --bracken_genus ${BASE}/RESULTS/BRACKEN/bracken.genus.mpa.report.txt \
    --bracken_species ${BASE}/RESULTS/BRACKEN/bracken.species.mpa.report.txt \
    --output ${BASE}/RESULTS/CONSENSUS_TAXA/bracken.metaphlan.taxa.prop.pdf \
    --output_common_genus ${BASE}/RESULTS/CONSENSUS_TAXA/bracken.metaphlan.common.genus.mpa.report.txt \
    --output_common_species ${BASE}/RESULTS/CONSENSUS_TAXA/bracken.metaphlan.common.species.mpa.report.txt

echo "=========================================="
echo "Consensus taxa complete!"
echo "Finished: $(date)"
echo "=========================================="
