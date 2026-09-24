#!/usr/bin/env python3
"""Alignment validation, step 1: which species get a reference genome.

Input: the species table going downstream (consensus species table, or the Bracken species table
when consensus passed through). A species is listed if it has >= --min_support reads in at least
one sample; below that it can never validate, so it isn't worth a genome.

Species are matched to NCBI taxids through the Kraken database's own taxonomy (taxDB: taxid,
parent, name, rank), so names are exactly those Bracken reports. Names that don't resolve are
kept with taxid NA and become `not_in_reference`.

Outputs:
  validation_species.tsv   label, species, taxid, max_reads, n_samples_at_support
  ref_key.txt              a hash of everything that determines the reference, used as its
                           cache key (sorted taxids + build settings)
"""
import argparse, hashlib, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_common import label_of, read_species_table, spaced, species_from_clade


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species_table", required=True)
    ap.add_argument("--taxdb", required=True, help="taxDB from the KrakenUniq database")
    ap.add_argument("--min_support", type=int, default=25)
    ap.add_argument("--key_settings", default="", help="build settings folded into the cache key")
    ap.add_argument("--out", default="validation_species.tsv")
    ap.add_argument("--key_out", default="ref_key.txt")
    a = ap.parse_args(argv)

    samples, table = read_species_table(a.species_table)
    best = {}
    for clade, counts in table.items():
        # Not candidates for a bacterial/archaeal/eukaryotic reference (audit B09): viruses (consensus
        # drops them, but a pass-through or --skip_consensus table does not) and human, which is
        # already in the index as T2T under the Homo_sapiens label.
        if "Viruses" in clade or label_of(species_from_clade(clade)) == "Homo_sapiens":
            continue
        name = species_from_clade(clade)
        mx = max(counts.values()) if counts else 0
        n_ok = sum(v >= a.min_support for v in counts.values())
        prev = best.get(name, (0, 0))
        best[name] = (max(prev[0], mx), max(prev[1], n_ok))
    wanted = {n: v for n, v in best.items() if v[0] >= a.min_support}

    # Match on the normalised label: the table says 'Escherichia_coli', taxDB 'Escherichia coli'.
    by_label = {label_of(n): n for n in wanted}
    taxid, canonical, ambiguous = {}, {}, set()
    with open(a.taxdb) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[3] == "species" and label_of(f[2]) in by_label:
                n = by_label[label_of(f[2])]
                if n not in taxid:
                    taxid[n], canonical[n] = f[0], f[2]
                elif taxid[n] != f[0]:
                    ambiguous.add(n)   # two taxDB species share this label (audit B01)
    # An ambiguous label gets no taxid rather than whichever came first in taxDB; the builder then
    # falls back to the name, and the species is flagged below.
    for n in ambiguous:
        del taxid[n]

    rows = sorted(wanted.items(), key=lambda kv: -kv[1][0])
    with open(a.out, "w") as out:
        out.write("label\tspecies\ttaxid\tmax_reads\tn_samples_at_support\n")
        for name, (mx, n_ok) in rows:
            # species column in taxDB/RefSeq spelling (spaces), for the name fallback in the build
            out.write(f"{label_of(name)}\t{canonical.get(name, spaced(name))}\t{taxid.get(name, 'NA')}\t{mx:g}\t{n_ok}\n")

    if not rows:
        sys.exit(f"ERROR: no species has >= {a.min_support} reads in any sample of {a.species_table}; "
                 "there is nothing to validate. Lower --validation_min_support or turn validation off.")
    key_src = ",".join(sorted(taxid.get(n, "NA:" + n) for n in wanted)) + "|" + a.key_settings
    with open(a.key_out, "w") as out:
        out.write(hashlib.sha256(key_src.encode()).hexdigest()[:16] + "\n")
    unresolved = [n for n in wanted if n not in taxid]
    print(f"{len(rows)} species with >= {a.min_support} reads in >= 1 sample; "
          f"{len(unresolved)} not resolved to a taxid ({len(ambiguous)} ambiguous: one label, several taxDB taxids)")


if __name__ == "__main__":
    main()
