#!/usr/bin/env bash
# Executed integration test for CMPipeline (audit C22). Unlike -preview and the unit tests, this
# RUNS the real main.nf, locally, on synthetic data only (random bases; see make_fixtures.py):
#
#   run A  --start_from beginning: BAM + CRAM (with a synthetic reference) + paired FASTQ x2 +
#          a single-end thin sample -> extraction, fastp, FastQC, host mapping against two tiny
#          synthetic minimap2 indexes, stub KrakenUniq + REAL Bracken and the MetaPhlAn gate,
#          consensus (Bracken pass-through), decontam on a 2-batch table, QC summary, MultiQC,
#          cohort report.
#   run B  --start_from decontam on a synthetic 12-sample, 2-batch table; batch correction too
#          when the batch-correction env exists, expecting the loud "Insufficient samples" stop.
#
# Usage: tests/integration/run_integration.sh <fixture_dir>
#   NEXTFLOW=<nextflow>          default: nextflow on PATH
#   CMP_IT_BATCH_CORR_ENV=<dir>  batch-correction env prefix; default --batch_corr_env from main.nf
#   CMP_IT_SKIP_RUNS=1           only regenerate fixtures and re-check existing outputs
# Prints PASS/FAIL per assertion and counts only. Exit status = number of failed assertions.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HERE="$REPO/tests/integration"
FX_ARG="${1:?usage: run_integration.sh <fixture_dir>}"
mkdir -p "$FX_ARG"
FX="$(cd "$FX_ARG" && pwd)"
NEXTFLOW="${NEXTFLOW:-nextflow}"
# The pipeline's own samtools and minimap2 envs build the fixtures (same versions the run uses).
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$REPO/.conda_cache/env-fcc8f6adcd9d7d83ce3eda3244306cdd/bin}"
MINIMAP2_BIN="${MINIMAP2_BIN:-$REPO/.conda_cache/env-8a53a5d3b3e246627eb5ffa09b5e2bcd/bin}"
BATCH_ENV="${CMP_IT_BATCH_CORR_ENV:-/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/envs/batch_correction_20260923}"

# Only ever wipe a directory this script created (marker file), never an arbitrary path.
MARK="$FX/.cmp_integration_fixture"
if [[ -n "$(ls -A "$FX" 2>/dev/null)" && ! -e "$MARK" ]]; then
    echo "FAIL  fixture dir is not empty and not a CMPipeline integration fixture dir: $FX" >&2
    exit 1
fi
touch "$MARK"

F="$FX/fixtures"
if [[ "${CMP_IT_SKIP_RUNS:-0}" != 1 ]]; then
    rm -rf "$F" "$FX/runA" "$FX/runB" "$FX/runC" "$FX/runD"
    mkdir -p "$F" "$FX/runA" "$FX/runB" "$FX/runC" "$FX/runD"
    if ! PATH="$SAMTOOLS_BIN:$MINIMAP2_BIN:$PATH" python3 "$HERE/make_fixtures.py" --outdir "$F" > "$FX/make_fixtures.log" 2>&1; then
        echo "FAIL  fixture generation (see $FX/make_fixtures.log)"
        exit 1
    fi
    echo "PASS  fixtures generated"

    common=(-c "$HERE/integration.config" -ansi-log false
            --control_role surrogate_biological --positive_values Tumor --control_values Blood
            --taxonomic_rank genus --run_decontam true)

    echo "INFO  run A (full pipeline) ..."
    ( cd "$FX/runA" && "$NEXTFLOW" run "$REPO/main.nf" "${common[@]}" \
        --sample "$F/samplesheet.csv" --outdir "$FX/runA/results" \
        --hg38_db "$F/db/hg38.mmi" --t2t_phix_db "$F/db/t2t_phix.mmi" \
        --kraken_db "$F/db/kraken" --metaphlan_db "$F/db/metaphlan" \
        --metadata_file "$F/metadata_runA.tsv" --min_total_per_batch 2 \
        --run_batch_correction false ) > "$FX/runA/nextflow.stdout" 2>&1
    echo $? > "$FX/runA/exit_status"
    echo "INFO  run A exit $(cat "$FX/runA/exit_status")"

    if [[ -d "$BATCH_ENV" ]]; then BC=on; else BC=off; fi
    echo "INFO  run B (--start_from decontam, batch correction $BC) ..."
    [[ $BC == off ]] && echo "SKIP  batch correction leg: env not found at $BATCH_ENV (set CMP_IT_BATCH_CORR_ENV)"
    ( cd "$FX/runB" && "$NEXTFLOW" run "$REPO/main.nf" "${common[@]}" \
        --start_from decontam --outdir "$FX/runB/results" \
        --consensus_otu_table "$F/decontam/genus_table.tsv" \
        --decontam_library_sizes "$F/decontam/library_sizes.tsv" \
        --metadata_file "$F/decontam/metadata.tsv" \
        --run_batch_correction "$([[ $BC == on ]] && echo true || echo false)" \
        --batch_corr_env "$BATCH_ENV" --batch_corr_covariates age_diag,sex,bmi ) > "$FX/runB/nextflow.stdout" 2>&1
    echo $? > "$FX/runB/exit_status"
    echo "$BC" > "$FX/runB/batch_correction"
    echo "INFO  run B exit $(cat "$FX/runB/exit_status")"

    # run C (audit C22): 28 samples, 2 batches, real batch shift -> batch correction must CORRECT.
    if [[ $BC == on ]]; then
        echo "INFO  run C (--start_from decontam + batch correction on 28 samples) ..."
        ( cd "$FX/runC" && "$NEXTFLOW" run "$REPO/main.nf" "${common[@]}" \
            --start_from decontam --outdir "$FX/runC/results" \
            --consensus_otu_table "$F/batch/genus_table.tsv" \
            --decontam_library_sizes "$F/batch/library_sizes.tsv" \
            --metadata_file "$F/batch/metadata.tsv" --run_batch_correction true \
            --batch_corr_env "$BATCH_ENV" --batch_corr_covariates age_diag,sex,bmi ) > "$FX/runC/nextflow.stdout" 2>&1
        echo $? > "$FX/runC/exit_status"
        echo "INFO  run C exit $(cat "$FX/runC/exit_status")"
    fi
    # run D (audit C04): a covariate missing from the metadata must stop at LAUNCH, before any task.
    echo "INFO  run D (missing covariate, expect a launch-time stop) ..."
    ( cd "$FX/runD" && "$NEXTFLOW" run "$REPO/main.nf" "${common[@]}" \
        --start_from decontam --outdir "$FX/runD/results" \
        --consensus_otu_table "$F/batch/genus_table.tsv" \
        --decontam_library_sizes "$F/batch/library_sizes.tsv" \
        --metadata_file "$F/batch/metadata_nocov.tsv" --run_batch_correction true \
        --batch_corr_env "$BATCH_ENV" --batch_corr_covariates age_diag,sex,bmi ) > "$FX/runD/nextflow.stdout" 2>&1
    echo $? > "$FX/runD/exit_status"
fi

python3 "$HERE/check_outputs.py" --fixtures "$F" --run-a "$FX/runA" --run-b "$FX/runB" \
    --run-c "$FX/runC" --run-d "$FX/runD" \
    --batch-correction "$(cat "$FX/runB/batch_correction" 2>/dev/null || echo off)"
