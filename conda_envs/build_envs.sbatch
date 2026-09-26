#!/usr/bin/env bash
#SBATCH --job-name=cmp-build-envs
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=build_envs-%j.log
#
# Builds CMPipeline's exact, locked conda environments (audit D01/D02/D03/D04, 2026-09-23).
# No data is touched. Run from anywhere:   sbatch conda_envs/build_envs.sbatch
#
#  1. envs whose yml did NOT change are locked from the build today's runs used (no rebuild)
#  2. metaphlan4 / humann3 are built fresh from their yml, then locked
#  3. (batch correction: see the note at step 3 below; Nextflow builds it)
# Lockfiles: conda_envs/lock/<env>.linux-64.txt  (`conda list --explicit --md5`)
set -euo pipefail
REPO=/tscc/lustre/restricted/alexandrov-ddn/users/shared/OncoBiomeProfiler-paper/workflow/CMPipeline
CACHE=$REPO/.conda_cache
LOCK=$REPO/conda_envs/lock
STAGE=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/cmp_env_build
source /tscc/projects/ps-lalexandrov/amabbasi/anaconda3/etc/profile.d/conda.sh
export CONDA_CHANNEL_PRIORITY=strict
# Never let ~/.local or a shell PYTHONPATH leak into the import checks (audit D13): job 12331453
# failed because the HUMAnN env's metaphlan imported a broken ~/.local pandas.
export PYTHONNOUSERSITE=1
unset PYTHONPATH
mkdir -p "$LOCK" "$STAGE"

lock() {  # lock <env prefix> <name>
    conda list --explicit --md5 -p "$1" > "$LOCK/$2.linux-64.txt"
    if grep -q "repo.anaconda.com" "$LOCK/$2.linux-64.txt"; then
        echo "ERROR: $2 contains packages from the 'defaults' channel" >&2; exit 1
    fi
    echo "locked $2: $(grep -c '^https' "$LOCK/$2.linux-64.txt") packages"
}

echo "== 1. unchanged envs: lock today's builds"
for pair in \
    env-fcc8f6adcd9d7d83ce3eda3244306cdd:samtools_env  env-ce688a4c172ecb243d52aba28db7fbf7:fastp_env \
    env-7a5eb82bdd08630942ad055fbe747db1:fastqc_env    env-8a53a5d3b3e246627eb5ffa09b5e2bcd:minimap2_env \
    env-f4ee76d686275fda7040d25d4e55387a:krakenUniq_bracken_env env-7fd621308a973b851c83b8a2c3560420:consensus_taxa_env \
    env-23bef169e8cc02b68c684eb859cf0da1:decontam_env  env-d9f413e80ec57f504ceefbb079753d2f:multiqc_env \
    env-a29f63f80ea4c50487fcf594750fc8da:validation_env; do
    lock "$CACHE/${pair%%:*}" "${pair##*:}"
done

echo "== 2. changed envs: build fresh, then lock"
# An env that already built AND was locked from this same yml is reused (REBUILD=1 forces a rebuild).
for name in metaphlan4_env humann3_env; do
    if [[ "${REBUILD:-0}" != 1 && -x "$STAGE/$name/bin/python" && -s "$LOCK/$name.linux-64.txt" \
          && "$LOCK/$name.linux-64.txt" -nt "$REPO/conda_envs/$name.yml" ]]; then
        echo "reusing $STAGE/$name (locked, yml unchanged)"
    else
        rm -rf "$STAGE/$name"
        conda env create -q -p "$STAGE/$name" -f "$REPO/conda_envs/$name.yml"
    fi
    lock "$STAGE/$name" "$name"
done
# Real import checks: a noarch package can land in the wrong site-packages (it did: humann
# under python3.12 in a 3.11 env) and the env still "solves".
"$STAGE/metaphlan4_env/bin/python" -c "import metaphlan; print('metaphlan import OK')"
"$STAGE/metaphlan4_env/bin/metaphlan" --version
"$STAGE/humann3_env/bin/python" -c "import humann, metaphlan; print('humann + metaphlan import OK')"
"$STAGE/humann3_env/bin/humann" --version
"$STAGE/humann3_env/bin/metaphlan" --version

# == 3. batch correction: no longer built here. Nextflow builds it from
# lock/batch_correction_env.linux-64.txt, whose two @LOCAL_CHANNEL@ lines are the r-cqrreg and
# r-conqur packages in conda_envs/local_channel/ (built from conda_envs/recipes/ by
# scripts/build_local_conda_packages.sh). If that lockfile is ever regenerated, keep those two
# lines, with the @LOCAL_CHANNEL@ prefix, at the end of the file.
echo "== done  $(date)"
