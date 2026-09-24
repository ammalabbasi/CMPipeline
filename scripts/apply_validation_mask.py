#!/usr/bin/env python3
"""Alignment validation, step 4: mask the Bracken species table and roll it up to genus.

Per sample, a species' Bracken count is set to 0 ONLY where its validation status is `failed`
(reads present, but they don't come from that genome). `validated`, `insufficient` and
`not_in_reference` keep their counts; decontam's own >= 5-read / >= 5%-prevalence filters handle
the weak ones. Group members (e.g. E. coli + Shigella) share their group's status.

Outputs:
  validated.species.tsv       masked species table (same layout as the input, clade_name first)
  validated.genus.tsv         genus roll-up of the MASKED species counts
  validation_matrix.tsv       species x sample status
  validation_stats.all.tsv    every sample's stats, concatenated
  validation_mask_report.tsv  per species: samples validated / failed / insufficient, reads removed
"""
import argparse, collections, glob, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_common import genus_clade, label_of, read_groups, read_species_table, species_from_clade


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species_table", required=True)
    ap.add_argument("--stats", nargs="+", required=True, help="<sample>.validation_stats.tsv files")
    ap.add_argument("--groups", default=None)
    ap.add_argument("--out_dir", default=".")
    a = ap.parse_args(argv)

    groups = read_groups(a.groups)
    samples, table = read_species_table(a.species_table)

    status = {}          # (sample, taxon) -> status
    header = None
    # A sample is "seen" if its stats FILE exists: an all-zero sample legitimately has a
    # header-only file (no taxon to score), and must not abort the cohort (audit B02).
    seen_samples = {os.path.basename(p)[:-len(".validation_stats.tsv")] for p in a.stats}
    with open(os.path.join(a.out_dir, "validation_stats.all.tsv"), "w") as allout:
        for path in sorted(a.stats):
            with open(path) as fh:
                h = fh.readline()
                if header is None:
                    header = h; allout.write(h)
                cols = h.rstrip("\n").split("\t")
                for line in fh:
                    allout.write(line)
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    status[(r["sample"], r["taxon"])] = r["status"]
    missing = [s for s in samples if s not in seen_samples]
    if missing:
        sys.exit(f"ERROR: no validation stats for {len(missing)} sample(s) of the species table, e.g. {missing[:3]}")

    masked = {}
    report = collections.defaultdict(lambda: collections.Counter())
    for clade, counts in table.items():
        tax = groups.get(label_of(species_from_clade(clade)), label_of(species_from_clade(clade)))
        row = {}
        for s in samples:
            # a zero count is "absent" whatever its group's status (audit B06): it was never scored
            st = "absent" if counts[s] == 0 else status.get((s, tax), "unscored")
            report[clade][st] += 1
            if st in ("failed", "unsupported"):   # unsupported: decision 29
                report[clade]["reads_removed"] += counts[s]
                row[s] = 0.0
            else:
                row[s] = counts[s]
        masked[clade] = row

    def fmt(v):
        # Counts are written exactly: `:g` rounded anything >= 1e6 to 6 significant digits (audit A15).
        return str(int(v)) if float(v).is_integer() else repr(float(v))

    def write_table(path, rows):
        with open(path, "w") as out:
            out.write("clade_name\t" + "\t".join(samples) + "\n")
            for clade, r in rows.items():
                out.write(clade + "\t" + "\t".join(fmt(r[s]) for s in samples) + "\n")

    write_table(os.path.join(a.out_dir, "validated.species.tsv"), masked)
    genus = collections.defaultdict(lambda: collections.Counter())
    for clade, r in masked.items():
        for s in samples:
            genus[genus_clade(clade)][s] += r[s]
    write_table(os.path.join(a.out_dir, "validated.genus.tsv"), genus)

    with open(os.path.join(a.out_dir, "validation_matrix.tsv"), "w") as out:
        out.write("clade_name\t" + "\t".join(samples) + "\n")
        for clade, counts in table.items():
            tax = groups.get(label_of(species_from_clade(clade)), label_of(species_from_clade(clade)))
            out.write(clade + "\t" + "\t".join(
                "absent" if counts[s] == 0 else status.get((s, tax), "unscored") for s in samples) + "\n")

    keys = ("validated", "failed", "unsupported", "insufficient", "not_in_reference", "absent", "unscored")
    with open(os.path.join(a.out_dir, "validation_mask_report.tsv"), "w") as out:
        out.write("clade_name\t" + "\t".join(f"n_{k}" for k in keys) + "\treads_removed\n")
        for clade, c in report.items():
            out.write(clade + "\t" + "\t".join(str(c[k]) for k in keys) + f"\t{fmt(c['reads_removed'])}\n")
    n_failed = sum(c["failed"] + c["unsupported"] for c in report.values())
    print(f"masked {n_failed} sample x species cell(s); "
          f"{sum(c['reads_removed'] for c in report.values()):g} Bracken reads removed")


if __name__ == "__main__":
    main()
