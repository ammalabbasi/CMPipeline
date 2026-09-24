#!/usr/bin/env bash
#SBATCH --job-name=cmp-chocophlan
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=fetch_chocophlan-%j.log
#
# Downloads ChocoPhlAn v201901_v31 (public, 16.5 GB), the nucleotide database HUMAnN 3.9 needs with
# MetaPhlAn 4 profiles (audit D05; Ammal 2026-09-23). The existing v296_201901b (MetaPhlAn 3 era) is
# left untouched. UniRef90 v201901b, already on disk, is the matching protein database.
set -euo pipefail
DB=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/databases/humann3
URL=http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
OUT=$DB/chocophlan_v201901_v31
mkdir -p "$OUT" && cd "$OUT"
[[ -s full_chocophlan.v201901_v31.tar.gz ]] || curl -fL --retry 5 -C - -o full_chocophlan.v201901_v31.tar.gz "$URL"
EXPECT=16502062909
GOT=$(stat -c %s full_chocophlan.v201901_v31.tar.gz)
[[ "$GOT" == "$EXPECT" ]] || { echo "ERROR: size $GOT != $EXPECT" >&2; exit 1; }
gzip -t full_chocophlan.v201901_v31.tar.gz
tar -xzf full_chocophlan.v201901_v31.tar.gz
md5sum full_chocophlan.v201901_v31.tar.gz > full_chocophlan.v201901_v31.tar.gz.md5
echo "files: $(ls | wc -l)"; du -sh .
echo "done $(date)"
