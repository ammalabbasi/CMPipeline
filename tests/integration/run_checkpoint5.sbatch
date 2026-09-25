#!/bin/bash -l
#SBATCH --job-name=cmp_checkpoint5
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/tscc/nfs/home/kohjy2000/decontam_contract_runs/checkpoint5_%j.sbatch.out
set -euo pipefail
REPO=/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/01.pipeline/CMPipeline/CMPipeline_github_260924_audit_fixes_kjy_v3
RUN=/tscc/nfs/home/kohjy2000/decontam_contract_runs/checkpoint5_v3_${SLURM_JOB_ID}
RUNNER=/tscc/nfs/home/kohjy2000/.cache/OncoBiomeProfiler/runner_f037c56dfc7faa9e
BATCH=/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/envs/batch_correction_20260923
mkdir "$RUN"
printf 'RUN_DIR=%s\n' "$RUN"
trap 'code=$?; printf "%s\n" "$code" > "$RUN/exit_status"' EXIT
export JAVA_HOME="$RUNNER"
export PATH="$RUNNER/bin:/tscc/nfs/home/kohjy2000/anaconda3/condabin:$PATH"
export NXF_VER=25.10.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$RUNNER/bin/nextflow" -version > "$RUN/nextflow_version.txt"
grep -q '25.10.4' "$RUN/nextflow_version.txt"
"$BATCH/bin/Rscript" "$REPO/tests/test_batch_contract.R" "$REPO" "$RUN/unit" > "$RUN/unit.log" 2>&1
python3 "$REPO/tests/integration/make_batch_fixture.py" "$RUN/fixtures"
cat > "$RUN/small.config" <<'CFG'
executor.cpus = 2
executor.memory = '16.GB'
process.cpus = 2
process.memory = '8.GB'
process.time = '45m'
params.max_cpus = 2
params.max_memory = '16.GB'
CFG
mkdir "$RUN/positive"
cd "$RUN/positive"
"$RUNNER/bin/nextflow" -log nextflow.log run "$REPO/main.nf" -profile local -c "$RUN/small.config" \
  --start_from batch_correction --run_batch_correction true --run_decontam false \
  --decontam_otu_table "$RUN/fixtures/counts.tsv" --metadata_file "$RUN/fixtures/metadata.csv" \
  --sample_id_column patient --type_column sample_type --batch_column shipment_batch \
  --batch_corr_covariates sample_type --taxon_column clade_name --taxonomic_rank species \
  --outdir "$RUN/positive/results" --skip_cohort_report true > stdout_stderr.log 2>&1
"$BATCH/bin/Rscript" "$REPO/tests/integration/assert_batch_output.R" "$REPO" "$RUN/fixtures" \
  "$RUN/positive/results/05_BATCH_CORRECTION/batch_corr_run" > "$RUN/assertions.log" 2>&1
mkdir "$RUN/negative12"
cd "$RUN/negative12"
set +e
"$RUNNER/bin/nextflow" -log nextflow.log run "$REPO/main.nf" -profile local -c "$RUN/small.config" \
  --start_from batch_correction --run_batch_correction true --run_decontam false \
  --decontam_otu_table "$RUN/fixtures/counts12.tsv" --metadata_file "$RUN/fixtures/metadata.csv" \
  --sample_id_column patient --type_column sample_type --batch_column shipment_batch \
  --batch_corr_covariates sample_type --taxon_column clade_name --taxonomic_rank species \
  --outdir "$RUN/negative12/results" --skip_cohort_report true > stdout_stderr.log 2>&1
code=$?
set -e
test "$code" -ne 0
grep -q 'Insufficient samples' stdout_stderr.log
if find results -name ConQuR_tuned.tsv -print -quit | grep -q .; then
  printf 'ERROR: negative12 published a corrected table\n'; exit 1
fi
printf 'PASS: existing tuned 20-sample guard retained\n' >> "$RUN/assertions.log"
printf 'PASS\n' > "$RUN/checkpoint5_status.txt"
cat "$RUN/assertions.log"
