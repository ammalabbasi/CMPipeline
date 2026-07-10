#!/bin/bash


# Load Conda
source ~/.bashrc

conda activate microbiome
echo "conda env activated"

cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/CMPipeline/out/ref/
echo "cd out dir"

mkdir krakenUniq
mkdir metaphlan

cd krakenUniq


echo "Downloading Kraken"
# KrakenUniq database
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/kuniq_microbialdb_minus_kdb.20230808.tgz
echo "Done file1"

wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/database.kdb
echo "Done file2"

cd ../metaphlan
# Metaphlan database
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.tar
echo "DONE"

