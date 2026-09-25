import csv
import sys
from pathlib import Path

root = Path(sys.argv[1]); root.mkdir(parents=True, exist_ok=True)
samples = [f'S{i:02d}' for i in range(1, 25)] + ['ZERO1', 'ZERO2', 'MISSING', 'UNKNOWN']
taxa = [f'k__Bacteria|g__Test{i}|s__Test{i}_fixture' for i in range(1, 9)]
metadata = []
for i, sample in enumerate(samples):
    batch = 'batch1' if i < 12 else 'batch2'
    kind = 'Tumor' if i % 2 == 0 else 'Normal'
    if sample == 'ZERO1': batch = 'batch1'
    if sample == 'MISSING': kind = ''
    if sample == 'UNKNOWN': batch = 'Unknown'
    metadata.append([sample, kind, batch])
with (root / 'metadata.csv').open('w', newline='') as f:
    w = csv.writer(f, quoting=csv.QUOTE_ALL); w.writerow(['patient', 'sample_type', 'shipment_batch']); w.writerows(metadata)
with (root / 'counts.tsv').open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['clade_name', *samples])
    for j, taxon in enumerate(taxa):
        counts = [0 if s.startswith('ZERO') else (20 + j * 7 + (i * (j + 3)) % 23 + (90 if i >= 12 and j < 4 else 0) + (15 if i % 2 == 0 and j % 2 == 0 else 0)) for i, s in enumerate(samples)]
        w.writerow([taxon, *counts])
subset = [i for i in range(6)] + [i for i in range(12, 18)]
with (root / 'counts12.tsv').open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['clade_name', *[samples[i] for i in subset]])
    for j, taxon in enumerate(taxa):
        w.writerow([taxon, *[20 + j * 7 + (i * (j + 3)) % 23 for i in subset]])
