#!/usr/bin/env python3
"""Aggregate only surviving species after decontamination; never retest genus."""
import argparse
import re
from count_tables import read_counts, write_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--taxon_column", default="clade_name")
    a = ap.parse_args()
    samples, species = read_counts(a.input, "s", delimiter=",", taxon_column=a.taxon_column)
    genus = {}
    for taxon, counts in species.items():
        parent = re.sub(r"\|s__[^|]+$", "", taxon)
        if not re.search(r"(?:^|\|)g__[^|]+$", parent):
            raise ValueError(f"Species has no explicit genus parent: {taxon}")
        row = genus.setdefault(parent, {s: 0 for s in samples})
        for s in samples:
            row[s] += counts[s]
    write_counts(a.output, samples, genus, delimiter=",", taxon_column=a.taxon_column)


if __name__ == "__main__":
    main()
