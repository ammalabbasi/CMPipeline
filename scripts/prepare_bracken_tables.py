#!/usr/bin/env python3
"""Convert hierarchical Bracken MPA to exact-ID rank tables and non-overlapping totals."""
import argparse
import csv
from count_tables import read_counts, write_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out_genus", required=True)
    ap.add_argument("--out_species", required=True)
    ap.add_argument("--out_library", required=True)
    a = ap.parse_args()
    samples, genus = read_counts(a.genus, "g")
    species_samples, species = read_counts(a.species, "s")
    if set(samples) != set(species_samples):
        raise ValueError("Bracken genus/species sample IDs differ")
    # Preserve the existing consensus exclusion, consistently on both entry routes.
    genus = {t: r for t, r in genus.items() if "Viruses" not in t}
    species = {t: r for t, r in species.items() if "Viruses" not in t}
    write_counts(a.out_genus, samples, genus)
    write_counts(a.out_species, samples, species)
    with open(a.out_library, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "bracken_genus_total", "bracken_species_total"])
        for s in samples:
            writer.writerow([s, sum(r[s] for r in genus.values()), sum(r[s] for r in species.values())])


if __name__ == "__main__":
    main()
