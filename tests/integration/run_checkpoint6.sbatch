#!/bin/bash -l
#SBATCH --job-name=cmp_readiness_v5
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/tscc/nfs/home/kohjy2000/decontam_contract_runs/checkpoint6_%j.sbatch.out
set -euo pipefail
REPO=/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/01.pipeline/CMPipeline/CMPipeline_github_260924_audit_fixes_kjy_v5
RUN=/tscc/nfs/home/kohjy2000/decontam_contract_runs/checkpoint6_v5_${SLURM_JOB_ID}
RUNNER=/tscc/nfs/home/kohjy2000/.cache/OncoBiomeProfiler/runner_f037c56dfc7faa9e
FRESH=/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/01.pipeline/CMPipeline/batch_env_checkpoint6_${SLURM_JOB_ID}
mkdir "$RUN"
printf 'RUN_DIR=%s\nFRESH_ENV=%s\n' "$RUN" "$FRESH"
trap 'code=$?; printf "%s\n" "$code" > "$RUN/exit_status"' EXIT
export JAVA_HOME="$RUNNER"
export PATH="$RUNNER/bin:/tscc/nfs/home/kohjy2000/anaconda3/condabin:$PATH"
export NXF_VER=25.10.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export R_LIBS="$FRESH/lib/R/library" R_LIBS_USER="$FRESH/lib/R/library" R_LIBS_SITE="$FRESH/lib/R/library"
bash "$REPO/scripts/setup_batch_env.sh" --prefix "$FRESH" \
  --conda /tscc/nfs/home/kohjy2000/anaconda3/bin/conda \
  --source-dir /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/conqur_offline > "$RUN/bootstrap.log" 2>&1
set +e
bash "$REPO/scripts/setup_batch_env.sh" --prefix "$FRESH" > "$RUN/bootstrap_refusal.log" 2>&1
code=$?
set -e
test "$code" -ne 0; grep -q 'Refusing to overwrite existing prefix' "$RUN/bootstrap_refusal.log"
"$RUNNER/bin/nextflow" -version > "$RUN/nextflow_version.txt"
grep -q '25.10.4' "$RUN/nextflow_version.txt"
"$FRESH/bin/Rscript" --vanilla "$REPO/tests/test_batch_contract.R" "$REPO" "$RUN/unit" > "$RUN/unit.log" 2>&1
"$FRESH/bin/Rscript" --vanilla "$REPO/tests/test_batch_reporting.R" "$REPO" "$RUN/reporting_unit" > "$RUN/reporting_unit.log" 2>&1
python3 "$REPO/tests/integration/make_batch_fixture.py" "$RUN/fixtures"
cat > "$RUN/small.config" <<'CFG'
executor.cpus = 2
executor.memory = '16.GB'
process { withName: BatchCorrection { cpus = 2; memory = '8.GB'; time = '45m' } }
CFG
for LABEL in positiveA positiveB negative12; do
  mkdir "$RUN/$LABEL"; cd "$RUN/$LABEL"
  INPUT="$RUN/fixtures/counts.tsv"
  if [[ "$LABEL" == negative12 ]]; then INPUT="$RUN/fixtures/counts12.tsv"; fi
  set +e
  "$RUNNER/bin/nextflow" -log nextflow.log run "$REPO/main.nf" -profile local -c "$RUN/small.config" \
    --start_from batch_correction --run_batch_correction true --run_decontam false \
    --decontam_otu_table "$INPUT" --metadata_file "$RUN/fixtures/metadata.csv" \
    --sample_id_column patient --type_column sample_type --batch_column shipment_batch \
    --batch_corr_covariates sample_type --taxon_column clade_name --taxonomic_rank species \
    --batch_corr_seed 314159 --batch_corr_env "$FRESH" \
    --outdir "$RUN/$LABEL/results" --skip_cohort_report true > stdout_stderr.log 2>&1
  code=$?
  set -e
  if [[ "$LABEL" == negative12 ]]; then
    test "$code" -ne 0; grep -q 'Insufficient samples' stdout_stderr.log
    if find results -name ConQuR_tuned.tsv -print -quit | grep -q .; then exit 1; fi
    printf 'PASS: tuned 20-sample guard retained\n' >> "$RUN/assertions.log"
  else
    test "$code" -eq 0
    "$FRESH/bin/Rscript" --vanilla "$REPO/tests/integration/assert_batch_output.R" "$REPO" "$RUN/fixtures" \
      "$RUN/$LABEL/results/05_BATCH_CORRECTION/batch_corr_run_species" >> "$RUN/assertions.log" 2>&1
  fi
done
python3 "$REPO/tests/integration/compare_batch_runs.py" \
  "$RUN/positiveA/results/05_BATCH_CORRECTION/batch_corr_run_species" \
  "$RUN/positiveB/results/05_BATCH_CORRECTION/batch_corr_run_species" \
  "$RUN/reproducibility_hashes.tsv" >> "$RUN/assertions.log" 2>&1
printf 'PASS\n' > "$RUN/checkpoint6_status.txt"
cat "$RUN/assertions.log"
