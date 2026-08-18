#!/bin/bash


# Load Conda
source ~/.bashrc

conda activate microbiome
echo "conda env activated"

cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/out/
echo "cd ref dir"

bash ../scripts/download_references.sh

echo "DONE"
