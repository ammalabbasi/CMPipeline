#!/usr/bin/env python3
"""Aggregate only surviving species after decontamination; never retest genus.

A species whose lineage has no genus directly above it (e.g. `d__Bacteria|s__unidentified_...`,
`s__synthetic_construct`) stays in the species table but cannot be summed into a genus. It is left
out of the genus table and listed, with its counts, in --excluded.
"""
import argparse
import csv
import re
from count_tables import read_counts, write_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--excluded", required=True, help="TSV of species left out of the genus table")
    ap.add_argument("--taxon_column", default="clade_name")
    a = ap.parse_args()
    samples, species = read_counts(a.input, "s", delimiter=",", taxon_column=a.taxon_column)
    genus, excluded = {}, []
    for taxon, counts in species.items():
        parent = re.sub(r"\|s__[^|]+$", "", taxon)
        if not re.search(r"(?:^|\|)g__[^|]+$", parent):
            reason = "genus is not the direct parent" if "g__" in parent else "no genus in lineage"
            excluded.append((taxon, reason, sum(counts.values()), sum(1 for s in samples if counts[s] > 0)))
            continue
        row = genus.setdefault(parent, {s: 0 for s in samples})
        for s in samples:
            row[s] += counts[s]
    if not genus:
        raise ValueError("No surviving species has a genus parent; the genus table would be empty")
    write_counts(a.output, samples, genus, delimiter=",", taxon_column=a.taxon_column)
    with open(a.excluded, "w", newline="") as fh:
        out = csv.writer(fh, delimiter="\t", lineterminator="\n")
        out.writerow(["taxon", "reason", "total_count", "samples_present"])
        out.writerows(excluded)
    print(f"Genus roll-up: {len(genus)} genera from {len(species) - len(excluded)} species; "
          f"{len(excluded)} species without a genus parent left out (see {a.excluded})")


if __name__ == "__main__":
    main()
