#!/bin/bash


# Load Conda
source ~/.bashrc

# Utility script, not part of the pipeline run (audit D09). Set TOOLS_ENV to an env that has the
# tools below; `microbiome` is the historical name and does not exist for every user.
conda activate "${TOOLS_ENV:-microbiome}" || { echo "ERROR: conda env '${TOOLS_ENV:-microbiome}' not found; rerun with TOOLS_ENV=<env with the tools this script needs>" >&2; exit 1; }
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

