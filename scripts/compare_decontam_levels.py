#!/usr/bin/env python3
"""Compare decontam's species-level and genus-level calls on the validated tables.

For every contaminant genus: which of its species were called at species level. For every
contaminant species: whether its genus was called. Agreement is the robust claim; a species call
whose genus was not called (or the reverse) is where the level matters.

Output: decontam_species_vs_genus.tsv
"""
import argparse, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_common import genus_clade


def contaminants(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    flag = next(c for c in ("final_contaminant", "contaminant") if c in rows[0]) if rows else None
    return {r["taxon_id"] for r in rows if flag and r[flag].upper() == "TRUE"}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species_final", required=True, help="<prefix>.contaminants_final.tsv, species run")
    ap.add_argument("--genus_final", required=True, help="<prefix>.contaminants_final.tsv, genus run")
    ap.add_argument("--out", default="decontam_species_vs_genus.tsv")
    a = ap.parse_args(argv)
    sp, ge = contaminants(a.species_final), contaminants(a.genus_final)
    with open(a.out, "w") as out:
        out.write("taxon\tlevel\tcalled_at_this_level\tcounterpart_called\tagreement\n")
        for s in sorted(sp):
            g = genus_clade(s)
            out.write(f"{s}\tspecies\tTRUE\t{str(g in ge).upper()}\t{'agree' if g in ge else 'species_only'}\n")
        for g in sorted(ge):
            sps = [s for s in sp if genus_clade(s) == g]
            out.write(f"{g}\tgenus\tTRUE\t{str(bool(sps)).upper()}\t{'agree' if sps else 'genus_only'}\n")
    print(f"species contaminants: {len(sp)}; genus contaminants: {len(ge)}")


if __name__ == "__main__":
    main()
