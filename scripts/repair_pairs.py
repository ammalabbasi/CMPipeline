#!/usr/bin/env python3
"""
Re-synchronize paired-end FASTQ files after independent host depletion.

Our host depletion pipeline processes R1 and R2 independently through minimap2,
which can result in different read counts. This script re-pairs them by keeping
only reads present in both files (matched by read name).

Usage:
    python scripts/repair_pairs.py \
        --r1 RESULTS/MAPPED_READS/sample.R1.fastq.gz \
        --r2 RESULTS/MAPPED_READS/sample.R2.fastq.gz \
        --out-dir RESULTS/MAPPED_READS_PAIRED/

Outputs:
    {out_dir}/{sample}.R1.paired.fastq.gz
    {out_dir}/{sample}.R2.paired.fastq.gz
"""

import argparse
import gzip
import os
import sys
from pathlib import Path


def parse_read_name(header: str) -> str:
    """Extract read name from FASTQ header, stripping /1 /2 and comments."""
    name = header.split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def index_fastq(filepath: str) -> dict:
    """Read a FASTQ file and return {read_name: (header, seq, plus, qual)}."""
    reads = {}
    opener = gzip.open if filepath.endswith(".gz") else open
    with opener(filepath, "rt") as f:
        while True:
            header = f.readline().rstrip("\n")
            if not header:
                break
            seq = f.readline().rstrip("\n")
            plus = f.readline().rstrip("\n")
            qual = f.readline().rstrip("\n")
            name = parse_read_name(header)
            reads[name] = (header, seq, plus, qual)
    return reads


def write_paired(reads: dict, names: list, outpath: str):
    """Write reads in order of names list to gzipped FASTQ."""
    with gzip.open(outpath, "wt", compresslevel=4) as f:
        for name in names:
            header, seq, plus, qual = reads[name]
            f.write(f"{header}\n{seq}\n{plus}\n{qual}\n")


def main():
    parser = argparse.ArgumentParser(description="Re-pair FASTQ files")
    parser.add_argument("--r1", required=True, help="R1 FASTQ file (gzipped)")
    parser.add_argument("--r2", required=True, help="R2 FASTQ file (gzipped)")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Derive output names
    r1_name = Path(args.r1).name.replace(".fastq.gz", ".paired.fastq.gz")
    r2_name = Path(args.r2).name.replace(".fastq.gz", ".paired.fastq.gz")
    out_r1 = os.path.join(args.out_dir, r1_name)
    out_r2 = os.path.join(args.out_dir, r2_name)

    # Skip if outputs already exist
    if os.path.exists(out_r1) and os.path.exists(out_r2):
        print(f"SKIP: outputs already exist for {Path(args.r1).name}")
        return

    print(f"Reading R1: {args.r1}")
    r1_reads = index_fastq(args.r1)
    print(f"  R1 reads: {len(r1_reads):,}")

    print(f"Reading R2: {args.r2}")
    r2_reads = index_fastq(args.r2)
    print(f"  R2 reads: {len(r2_reads):,}")

    # Find shared read names, preserve R1 order
    shared = set(r1_reads.keys()) & set(r2_reads.keys())
    # Maintain original order from R1
    ordered_names = [name for name in r1_reads if name in shared]

    dropped_r1 = len(r1_reads) - len(shared)
    dropped_r2 = len(r2_reads) - len(shared)
    print(f"  Shared pairs: {len(shared):,}")
    print(f"  Dropped from R1: {dropped_r1:,}")
    print(f"  Dropped from R2: {dropped_r2:,}")

    if len(shared) == 0:
        print("  WARNING: No shared reads found! Skipping.")
        return

    print(f"Writing paired R1: {out_r1}")
    write_paired(r1_reads, ordered_names, out_r1)

    print(f"Writing paired R2: {out_r2}")
    write_paired(r2_reads, ordered_names, out_r2)

    print(f"  Done. {len(shared):,} paired reads written.")


if __name__ == "__main__":
    main()
