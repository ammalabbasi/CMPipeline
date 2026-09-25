"""Exact deterministic table comparison, excluding elapsed times and run metadata."""
from pathlib import Path
import csv
import hashlib
import math
import sys

a, b, out = map(Path, sys.argv[1:])
files = ['corrected/ConQuR_tuned.tsv', 'batch_correction_quality.tsv',
         'batch_correction_methods.tsv', 'batch_correction_run_settings.tsv',
         'batch_correction_sample_status.tsv', 'batch_correction_taxon_status.tsv',
         'permanova_summary.tsv']
# Output directory name == the batch script's --prefix (<run>_<rank>), which names the normalized files.
if a.name != b.name: raise SystemExit(f'FAIL: compare the same run/rank directories: {a.name} vs {b.name}')
files += [f'normalized/{a.name}_{source}_{method}.tsv'
          for source in ['raw','ConQuR_tuned'] for method in ['clr_pseudo','clr_czm','rclr']]
rows=[]
for name in files:
    x=(a/name).read_bytes(); y=(b/name).read_bytes()
    if x != y: raise SystemExit(f'FAIL: independent table differs: {name}')
    rows.append([name, hashlib.sha256(x).hexdigest(), hashlib.sha256(y).hexdigest(), 'identical'])
def records(path):
    with path.open() as f: return list(csv.DictReader(f, delimiter='\t'))
settings={r['key']:r['value'] for r in records(a/'batch_correction_run_settings.tsv')}
assert settings['seed']=='314159' and settings['threads']=='2'
assert settings['rng_kind']=="L'Ecuyer-CMRG"
q=records(a/'batch_correction_quality.tsv')[0]
assert q['execution_status']=='completed'
assert q['taxon_fit_diagnostics']=='not_exposed_by_pinned_dependency'
assert q['biological_signal_preservation']=='not_assessed'
if float(q['final_R2'])<0 or float(q['baseline_R2'])<=0:
    assert q['reported_reduction_percent']=='NA' and q['quality_status'].startswith('not_interpretable')
else:
    assert math.isfinite(float(q['reported_reduction_percent']))
perms=records(a/'permanova_summary.tsv')
assert all(r['N']=='24' and r['Evaluation'] in ('margin','terms_batch_last') for r in perms)
assert all(r['Distance']==('bray' if r['Transform']=='counts' else 'euclidean') for r in perms)
with out.open('w', newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['file','sha256_A','sha256_B','comparison']);w.writerows(rows)
print(f'PASS: {len(files)} independent output tables byte-identical, seed/threads recorded, quality labels and evaluation metadata verified')
