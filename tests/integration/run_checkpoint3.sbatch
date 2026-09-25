#!/bin/bash
#SBATCH --job-name=cmp_checkpoint3_v2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --account=ddp302
set -euo pipefail
: "${SLURM_JOB_ID:?Run this test on an allocated compute node}"
ROOT="/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/01.pipeline/CMPipeline/CMPipeline_github_260924_audit_fixes_kjy_v2"
OUT="/tscc/nfs/home/kohjy2000/decontam_contract_runs/checkpoint3_v2_${SLURM_JOB_ID}"
mkdir -p "$OUT"
exec > >(tee "$OUT/test.log") 2>&1
printf 'RUN_DIR=%s\nHOST=%s\n' "$OUT" "$(hostname)"
CACHE="/tscc/lustre/restricted/alexandrov-ddn/users/shared/OncoBiomeProfiler-paper/workflow/CMPipeline/.conda_cache"
R_ENV="$CACHE/decontam_env.linux-64-fd0954b62149383eaa4d7a8f17819ec8"
C_ENV=""
for e in "$CACHE"/consensus_taxa_env.linux-64-*; do
    if [ -x "$e/bin/python" ]; then C_ENV="$e"; break; fi
done
if [ -z "$C_ENV" ] || [ ! -x "$R_ENV/bin/Rscript" ]; then
    printf 'ERROR: required package environments unavailable; no bypass or environment build attempted\n' >&2
    exit 1
fi
export CMP_TEST_RSCRIPT="$R_ENV/bin/Rscript"
PY="$C_ENV/bin/python"
"$CMP_TEST_RSCRIPT" -e 'stopifnot(as.character(packageVersion("decontam")) == "1.22.0"); print(sessionInfo())'
printf 'DECONTAM_ENV=%s\nCONSENSUS_ENV=%s\n' "$R_ENV" "$C_ENV"
export CMP_TEST_ROOT="$ROOT"
"$PY" - <<'PY'
import os, sys, unittest
from pathlib import Path
root = Path(os.environ["CMP_TEST_ROOT"])
sys.path.insert(0, str(root / "tests"))
suite = unittest.TestSuite()
for name in ("test_decontam_zero_total", "test_checkpoint3_contracts"):
    suite.addTests(unittest.defaultTestLoader.loadTestsFromName(name))
result = unittest.TextTestRunner(verbosity=2).run(suite)
if result.skipped:
    print("ERROR: skipped tests do not count as validation:", result.skipped)
sys.exit(0 if result.wasSuccessful() and not result.skipped else 1)
PY
"$PY" "$ROOT/tests/integration/checkpoint3_fixture.py" create "$OUT/fixtures"
RUNNER="/tscc/nfs/home/kohjy2000/.cache/OncoBiomeProfiler/runner_f037c56dfc7faa9e"
export JAVA_HOME="$RUNNER" PATH="$RUNNER/bin:$PATH" NXF_VER=25.10.4
nextflow -version
mkdir -p "$OUT/module"
cd "$OUT/module"
nextflow -log "$OUT/module/nextflow.log" run "$ROOT/main.nf" -ansi-log false \
    -work-dir "$OUT/module/work" --outdir "$OUT/module/results" --tracedir "$OUT/module/pipeline_info" \
    --start_from decontam --run_decontam true --run_batch_correction false \
    --consensus_otu_table "$OUT/fixtures/otu.tsv" --metadata_file "$OUT/fixtures/metadata.csv" \
    --decontam_library_sizes "$OUT/fixtures/libraries.tsv" --taxonomic_rank species \
    --sample_id_column patient --type_column sample_type --batch_column shipment_batch \
    --control_role surrogate_biological --positive_values Tumor --control_values Normal \
    --decontam_env "$R_ENV" --consensus_taxa_env "$C_ENV" \
    --max_cpus 2 --max_memory 16.GB --skip_multiqc true --skip_cohort_report true
"$PY" "$ROOT/tests/integration/checkpoint3_fixture.py" check "$OUT/module/results"
mkdir -p "$OUT/entry_preview"
cd "$OUT/entry_preview"
nextflow -log "$OUT/entry_preview/nextflow.log" run "$ROOT/main.nf" -preview -ansi-log false \
    --sample "$OUT/fixtures/samples.csv" --skip_host_depletion true \
    --run_decontam false --skip_multiqc true --skip_cohort_report true \
    --kraken_db "$OUT/fixtures/kraken" --metaphlan_db "$OUT/fixtures/metaphlan" --metaphlan_index synthetic
printf 'PASS: checkpoint 3 bounded tests; upstream classification entry was previewed, not executed\n'
printf 'PASS\n' > "$OUT/checkpoint3_status.txt"
