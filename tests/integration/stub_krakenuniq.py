#!/usr/bin/env python3
"""Stand-in for `krakenuniq` in the executed integration test (audit C22). Test use only.

There is no toy KrakenUniq database, so this classifies the synthetic reads by their NAME,
which make_fixtures.py sets to mic_s<species taxid>_N or mic_g<genus taxid>_N (anything else
is unclassified), and writes what the real tool writes: the 9-column KrakenUniq report
(% reads taxReads kmers dup cov taxID rank taxName, names indented two spaces per level),
the per-read output, and the classified / unclassified FASTA. Real Bracken then runs on the
report. Accepts the arguments Modules/Bracken.nf passes. Standard library only.
"""

import argparse
import re
import sys
from pathlib import Path

NAME = re.compile(r"^mic_([sg])(\d+)_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--threads")
    ap.add_argument("--report-file", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--classified-out", required=True)
    ap.add_argument("--unclassified-out", required=True)
    ap.add_argument("reads", nargs="+")
    args = ap.parse_args()

    tax = {}
    for line in (Path(args.db) / "stub_taxonomy.tsv").read_text().splitlines():
        taxid, parent, rank, name = line.split("\t")
        tax[int(taxid)] = (int(parent), rank, name)
    children = {}
    for t, (parent, _, _) in tax.items():
        children.setdefault(parent, []).append(t)

    own = {}          # taxid -> reads assigned directly
    total = 0
    with open(args.output, "w") as out, open(args.classified_out, "w") as cls, \
            open(args.unclassified_out, "w") as uncls:
        for path in args.reads:
            with open(path) as fh:
                while True:
                    header = fh.readline()
                    if not header:
                        break
                    seq = fh.readline().rstrip("\n")
                    fh.readline()
                    fh.readline()
                    name = header[1:].split()[0]
                    total += 1
                    m = NAME.match(name)
                    taxid = int(m.group(2)) if m and int(m.group(2)) in tax else 0
                    if taxid:
                        own[taxid] = own.get(taxid, 0) + 1
                        cls.write(f">{name}\n{seq}\n")
                    else:
                        uncls.write(f">{name}\n{seq}\n")
                    out.write(f"{'C' if taxid else 'U'}\t{name}\t{taxid}\t{len(seq)}\t{taxid}:1\n")

    def clade(t):
        return own.get(t, 0) + sum(clade(c) for c in children.get(t, []))

    def pct(n):
        return f"{100.0 * n / total:.4f}" if total else "0.0000"

    rows = []
    classified = clade(1)
    rows.append((pct(total - classified), total - classified, total - classified, 0, 0, "NA", 0, "no rank", "unclassified"))

    def walk(t, depth):
        n = clade(t)
        if n == 0:
            return
        _, rank, name = tax[t]
        rows.append((pct(n), n, own.get(t, 0), n * 50, 1.0, 0.01, t, rank, "  " * depth + name))
        for c in sorted(children.get(t, [])):
            walk(c, depth + 1)

    walk(1, 0)
    with open(args.report_file, "w") as rep:
        rep.write("# KrakenUniq stub for the CMPipeline integration test (synthetic data)\n")
        rep.write("# CL: stub_krakenuniq.py\n\n")
        rep.write("%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n")
        for r in rows:
            rep.write("\t".join(map(str, r)) + "\n")
    print(f"stub krakenuniq: {total} reads, {classified} classified", file=sys.stderr)


if __name__ == "__main__":
    main()
