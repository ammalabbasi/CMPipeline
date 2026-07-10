#!/bin/bash

# Load Conda
source ~/.bashrc

# Activate nextflow environment
conda activate microbiome
echo "Nextflow microbiome environment activated"

# Change to pipeline directory
cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline_updated/CMPipeline_github_251124_update_kjy/
echo "Working directory: $(pwd)"

# Resume nextflow with email notification
echo "Resuming Nextflow pipeline..."
nextflow run main.nf \
    -resume \
    -N l1joseph@ucsd.edu

echo "Pipeline completed"
