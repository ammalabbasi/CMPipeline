"""Strict count-table boundaries shared by preprocessing and alignment validation."""
import csv
import re
from decimal import Decimal, InvalidOperation


def count_value(value):
    try:
        n = Decimal(value)
    except (InvalidOperation, ValueError):
        raise ValueError(f"Invalid count: {value!r}")
    if not n.is_finite() or n < 0 or n != n.to_integral_value():
        raise ValueError(f"Expected finite nonnegative integer count: {value!r}")
    return int(n)


def read_counts(path, rank, delimiter="\t", taxon_column="clade_name"):
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader, [])
        if len(header) < 2 or header[0] not in (taxon_column, "#Classification"):
            raise ValueError(f"{path}: expected taxonomy column and sample columns")
        raw = header[0] == "#Classification"
        samples = header[1:]
        if raw:
            suffix = re.compile(r"\.bracken\." + rank.upper() + r"(\.mpa)?\.krakenreport\.txt$")
            if any(not suffix.search(s) for s in samples):
                raise ValueError(f"{path}: unrecognized Bracken sample label")
            samples = [suffix.sub("", s) for s in samples]
        if any(not s or s != s.strip() for s in samples) or len(set(samples)) != len(samples):
            raise ValueError(f"{path}: empty, padded or duplicate sample IDs")
        rows, seen = {}, set()
        terminal = re.compile(r"(?:^|\|)" + rank + r"__[^|]+$")
        for fields in reader:
            if len(fields) != len(header):
                raise ValueError(f"{path}: wrong field count on row {reader.line_num}")
            taxon = fields[0]
            if not taxon or taxon != taxon.strip() or taxon in seen:
                raise ValueError(f"{path}: empty, padded or duplicate taxon: {taxon!r}")
            seen.add(taxon)
            counts = {s: count_value(v) for s, v in zip(samples, fields[1:])}
            if not terminal.search(taxon):
                if raw:
                    continue
                raise ValueError(f"{path}: taxon is not rank {rank}: {taxon}")
            rows[taxon] = counts
    return samples, rows


def write_counts(path, samples, rows, delimiter="\t", taxon_column="clade_name"):
    with open(path, "w", newline="") as fh:
        out = csv.writer(fh, delimiter=delimiter, lineterminator="\n")
        out.writerow([taxon_column] + samples)
        for taxon, values in rows.items():
            out.writerow([taxon] + [values[s] for s in samples])
