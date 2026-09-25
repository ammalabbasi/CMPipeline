#!/usr/bin/env python3
"""Synthetic species/metadata inputs and exact post-decontamination genus assertions."""
import csv
import gzip
import sys
from pathlib import Path

mode, directory = sys.argv[1:3]
root = Path(directory)
ids = [f"S{i}" for i in range(1, 13)]
positive = [i % 6 < 3 for i in range(12)]
if mode == "create":
    root.mkdir(parents=True, exist_ok=True)
    with open(root / "otu.tsv", "w") as f:
        f.write("clade_name\t" + "\t".join(ids) + "\n")
        for taxon, counts in [
            ("d__Bacteria|g__G|s__G_positive", [10 if p else 0 for p in positive]),
            ("d__Bacteria|g__G|s__G_control", [0 if p else 40 for p in positive]),
            ("d__Bacteria|g__Other|s__Other_shared", [10] * 12)
        ]:
            f.write(taxon + "\t" + "\t".join(map(str, counts)) + "\n")
    with open(root / "metadata.csv", "w", newline="") as f:
        out = csv.writer(f, quoting=csv.QUOTE_ALL)
        out.writerow(["patient", "sample_type", "shipment_batch"])
        for i, s in enumerate(ids):
            out.writerow([s, "Tumor" if positive[i] else "Normal", f"batch{i // 6 + 1}"])
    with open(root / "libraries.tsv", "w") as f:
        f.write("sample\tbracken_genus_total\tbracken_species_total\n")
        for s in ids:
            f.write(f"{s}\t5000\t5000\n")
    fq = root / "host_depleted.fastq.gz"
    with gzip.open(fq, "wt") as f:
        f.write("@synthetic/1\n" + "ACGT" * 25 + "\n+\n" + "I" * 100 + "\n")
    with open(root / "samples.csv", "w", newline="") as f:
        out = csv.writer(f)
        out.writerow(["patient", "fastq_1"])
        out.writerow(["S1", str(fq)])
    for name in ["kraken", "metaphlan"]:
        (root / name).mkdir(exist_ok=True)
    (root / "kraken/taxDB").write_text("synthetic preview only\n")
    (root / "kraken/database50mers.kmer_distrib").write_text("synthetic preview only\n")
    (root / "metaphlan/synthetic.pkl").write_text("synthetic preview only\n")
elif mode == "check":
    d = root / "04_DECONTAMINATION/decontam_run"
    with open(d / "decontam_run.decontam_pkg_decontaminated.csv") as f:
        species = list(csv.DictReader(f))
    with open(d / "decontam_run.decontam_pkg_decontaminated.genus.csv") as f:
        genus = list(csv.DictReader(f))
    assert set(r["clade_name"] for r in species) == {
        "d__Bacteria|g__G|s__G_positive", "d__Bacteria|g__Other|s__Other_shared"
    }, species
    assert set(r["clade_name"] for r in genus) == {"d__Bacteria|g__G", "d__Bacteria|g__Other"}, genus
    for row in genus:
        for i, s in enumerate(ids):
            wanted = (10 if positive[i] else 0) if row["clade_name"].endswith("|g__G") else 10
            assert int(row[s]) == wanted, (row, s, wanted)
    with open(d / "decontam_run.validated_metadata.tsv") as f:
        assert [r["patient"] for r in csv.DictReader(f, delimiter="\t")] == ids
    print("PASS: official species decontamination, quoted CSV and exact post-removal genus sums")
else:
    raise ValueError("Expected create or check")
