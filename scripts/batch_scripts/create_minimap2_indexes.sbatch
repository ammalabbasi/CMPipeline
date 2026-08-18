#!/bin/bash


# Load Conda
source ~/.bashrc

conda activate microbiome
echo "conda env activated"

cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/out/
echo "cd out dir"

bash ../scripts/create_minimap2_indexes.sh

echo "DONE"

