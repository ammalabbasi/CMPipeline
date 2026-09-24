#!/usr/bin/env python3
"""Compare CRAM @SQ sequence names, lengths, and MD5s with a FASTA dictionary."""

import argparse
import subprocess
import sys


def sq_records(command):
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    records = {}
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        if "SN" in fields:
            if fields["SN"] in records:
                raise ValueError(f"duplicate @SQ sequence name: {fields['SN']}")
            records[fields["SN"]] = (fields.get("LN"), fields.get("M5"))
    return records


def validate(alignment, reference):
    cram = sq_records(["samtools", "view", "-H", alignment])
    ref = sq_records(["samtools", "dict", reference])
    if not cram or not ref:
        raise ValueError("CRAM header or reference dictionary contains no @SQ records")

    errors = []
    for name, values in cram.items():
        if name not in ref:
            errors.append(f"missing reference contig {name}")
        elif values != ref[name]:
            errors.append(f"{name}: CRAM LN/M5={values}, reference LN/M5={ref[name]}")
    if errors:
        raise ValueError("CRAM/reference mismatch: " + "; ".join(errors[:10]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--reference", required=True)
    args = parser.parse_args()
    validate(args.alignment, args.reference)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
