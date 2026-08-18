#!/bin/bash


# Load Conda
source ~/.bashrc

conda activate microbiome
echo "conda env activated"

cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/out/ref/
echo "cd out dir"

echo "unpacking kraken"
tar -xzf krakenUniq/kuniq_microbialdb_minus_kdb.20230808.tgz -C krakenUniq/

echo "unpacking metaphlan"
tar -xf metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202307.tar -C metaphlan/
echo "DONE"

