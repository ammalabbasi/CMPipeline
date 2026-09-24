#!/usr/bin/env bash
set -euo pipefail

# Fast, dependency-light checks for Module 1/2 wiring and portable environments.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

python3 - <<'PY'
from pathlib import Path
import yaml
files = sorted(Path('conda_envs').glob('*.yml'))
for path in files:
    spec = yaml.safe_load(path.read_text())
    assert isinstance(spec, dict), path
    assert 'prefix' not in spec, path
print(f'Conda YAML checks passed ({len(files)} files)')
PY

python3 -m py_compile scripts/compute_consensus_taxa.py scripts/validate_cram_reference.py
python3 -m unittest discover -s tests -p 'test_*.py'
# Executed integration test on synthetic data (audit C22). Opt-in: minutes, and needs the conda envs.
if [[ -n "${CMP_INTEGRATION_FIXTURES:-}" ]]; then
    bash tests/integration/run_integration.sh "$CMP_INTEGRATION_FIXTURES"
fi
git diff --check

if command -v nextflow >/dev/null 2>&1; then
    nextflow config -profile local >/dev/null
    nextflow config -profile tscc >/dev/null
    echo 'Nextflow profile checks passed (local, tscc)'
else
    echo 'Nextflow not found; skipped profile checks' >&2
fi

echo 'CMPipeline smoke checks passed'
