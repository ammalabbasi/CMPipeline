#!/usr/bin/env python3
"""Input provenance from a BAM/CRAM header: which reference, which aligner chain (audit X13).

Ludmil E3 / pksProfiler 70cb629 (scripts/header_provenance.py), adapted. A results tree could
not say which human build each BAM/CRAM was aligned to, nor how. That matters for a microbiome
read-out from unmapped reads: TCGA (GDC GRCh38.d1.vd1: hs38d1 decoys, chrEBV and ~200 viral
contigs), Mutographs (1000G GRCh38 + decoys/HLA) and Hartwig (alt-masked GRCh38) differ in the
contig classes reads can be absorbed by, so they differ in which reads FLAG 4 extraction can see.

extractReads already writes the header once (`samtools view -H` -> input.header.sam) for the
low-score rescue decision; scripts/low_score_rescue.py decide --provenance calls parse() on that
same text, so the alignment header is read once per sample. Header lines only; no alignment record
is read.

One row per sample (PROVENANCE_FIELDS):

  input_ref_fingerprint     md5 over the ordered "SN:LN" lines of @SQ. Same reference, same
                            order -> same value; a different build or contig set -> different
  input_ref_m5_fingerprint  md5 over the sorted, de-duplicated @SQ M5 values; NA unless every @SQ
                            carries M5. Insensitive to contig renaming (chr1 vs 1)
  input_ref_contigs         number of @SQ lines
  input_ref_build           GRCh38 / GRCh37 / T2T-CHM13 from chr1's length; unknown; no_chr1
  input_ref_alt_contigs, input_ref_decoy_contigs, input_ref_hla_contigs, input_ref_viral_contigs
                            @SQ counts per contig class (CONTIG_CLASSES; names only). Viral
                            includes chrEBV
  input_ref_has_ebv         true / false: chrEBV (or NC_007605) is a contig of the reference
  input_pg_chain            the @PG records in PP order as PN:VN (ID when no PN), " > "-joined.
                            Program names and versions only: the CL is never written
  input_bwa_options         the bwa mem / bwa-mem2 mem options that change alignments, reduced to
                            option letters and NUMERIC values only (e.g. "-K 100000000 -M -T 0 -Y").
                            Never -R (read group), -H, -o, file names or any other text: those
                            carry sample barcodes and paths. NA without a bwa mem @PG
  input_unmapped_filter     none, or the upstream step(s) in the @PG chain that removed unmapped
                            reads before this file was written (FILTER_REASONS): samtools view
                            -F with bit 4, -f with bit 2, -q > 0, -L/region; sambamba view
                            -F "...unmapped..."/mapping_quality. A filtered input holds no (or
                            not all) unmapped reads, so an empty profile from it is not a
                            microbial negative (QC floor, build_qc_summary.py terminal_state)

The cohort aggregate (finalize_provenance.py, RUN_REPORT "INPUT ALIGNMENTS") is counts only.
Standard library only.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import re
import shlex
from pathlib import Path

NA = "NA"

PROVENANCE_FIELDS = [
    "sample", "input_ref_fingerprint", "input_ref_m5_fingerprint", "input_ref_contigs",
    "input_ref_build", "input_ref_alt_contigs", "input_ref_decoy_contigs", "input_ref_hla_contigs",
    "input_ref_viral_contigs", "input_ref_has_ebv", "input_pg_chain", "input_bwa_options",
    "input_unmapped_filter",
]

# chr1 length -> build. The GRCh38 analysis sets (GDC d1.vd1, 1000G hs38DH, alt-masked) all share it.
CHR1_BUILD = {248956422: "GRCh38", 249250621: "GRCh37", 248387328: "T2T-CHM13"}

# Contig classes, from the @SQ name only. First match wins.
CONTIG_CLASSES = (
    # EBV is counted as viral and also reported on its own (input_ref_has_ebv).
    ("viral", re.compile(r"^(chrEBV|NC_007605|EBV|phi ?X|NC_001422|HPV|HBV|HCV|HIV|HTLV|KSHV|HHV|CMV|"
                         r"MCV|MCPyV|SV40)", re.I)),
    ("hla", re.compile(r"^HLA-", re.I)),
    ("alt", re.compile(r"_alt$", re.I)),
    ("decoy", re.compile(r"_decoy$|^chrUn_JTFH|^hs37d5$", re.I)),
)
EBV = re.compile(r"^(chrEBV|NC_007605|EBV)$", re.I)

# --- bwa mem options worth recording (bwa 0.7.x and bwa-mem2 share them). Everything else in the
# CL is dropped. Values are kept only when they are numbers (or number lists, e.g. -O 6,6).
BWA_PROGRAMS = ("bwa", "bwa-mem2")
BWA_NUMERIC = set("kwdrycDWmABOELUTKIhzs")   # options taking a numeric (or numeric-list) value
BWA_FLAGS = set("SPpaCVYMj5qe")              # options without a value
BWA_SKIP_VALUE = set("RHotvxGF")             # options with a value that is never recorded
NUMERIC_VALUE = re.compile(r"^-?\d+(\.\d+)?([eE][-+]?\d+)?([,:]-?\d+(\.\d+)?)*$")

# input_unmapped_filter vocabulary. build_qc_summary.py and the docs use exactly these.
FILTER_REASONS = ("samtools_view_exclude_unmapped", "samtools_view_require_proper_pair",
                  "samtools_view_min_mapq", "samtools_view_region", "sambamba_view_filter")
SAMTOOLS_FLAG_NAMES = {
    "PAIRED": 0x1, "PROPER_PAIR": 0x2, "UNMAP": 0x4, "MUNMAP": 0x8, "REVERSE": 0x10,
    "MREVERSE": 0x20, "READ1": 0x40, "READ2": 0x80, "SECONDARY": 0x100, "QCFAIL": 0x200,
    "DUP": 0x400, "SUPPLEMENTARY": 0x800,
}


def tags(line: str) -> dict:
    out = {}
    for field in line.rstrip("\n").split("\t")[1:]:
        if len(field) > 3 and field[2] == ":":
            out.setdefault(field[:2], field[3:])
    return out


def md5(lines) -> str:
    return hashlib.md5("\n".join(lines).encode()).hexdigest()


def clean(value) -> str:
    """A header value safe for a TSV cell: no tabs/newlines/spaces, bounded length."""
    if not value:
        return NA
    return re.sub(r"[^A-Za-z0-9._+:-]", "_", str(value))[:48] or NA


def _program(token: str) -> str:
    return token.strip("'\"").rsplit("/", 1)[-1]


def _split(cl: str) -> list[str]:
    try:
        return shlex.split(cl)
    except ValueError:            # unbalanced quotes in a free-text CL
        return cl.split()


def contig_class(name: str) -> str:
    for label, pattern in CONTIG_CLASSES:
        if pattern.search(name):
            return label
    return "other"


def pg_chain(pgs: list[dict]) -> list[dict]:
    """@PG records ordered along their PP links, roots first; unreachable ones appended."""
    by_id = {p.get("ID"): p for p in pgs}
    children = collections.defaultdict(list)
    for p in pgs:
        children[p.get("PP")].append(p)
    roots = [p for p in pgs if p.get("PP") not in by_id] or pgs[:1]
    out, seen = [], set()
    for root in roots:
        stack = [root]
        while stack:
            node = stack.pop()
            if node.get("ID") in seen:
                continue
            seen.add(node.get("ID"))
            out.append(node)
            stack.extend(reversed(children.get(node.get("ID"), [])))
    out.extend(p for p in pgs if p.get("ID") not in seen)
    return out


def pg_label(pg: dict) -> str:
    """PN:VN; with no PN, the ID without the .N / -N / -HEX suffix samtools and Picard add."""
    name = pg.get("PN") or re.sub(r"([.-]\d+|-[0-9A-F]{6,})$", "", pg.get("ID", ""), flags=re.I)
    return f"{clean(name)}:{clean(pg.get('VN'))}"


MAX_CHAIN = 12


def chain_text(pgs: list[dict]) -> str:
    if not pgs:
        return NA
    labels = [pg_label(p) for p in pg_chain(pgs)]
    if len(labels) > MAX_CHAIN:
        labels = labels[:MAX_CHAIN] + [f"+{len(labels) - MAX_CHAIN}_more"]
    return " > ".join(labels)


def bwa_mem_options(cl: str):
    """The reduced bwa mem option string of one CL, or None if it does not run bwa mem."""
    tokens = _split(cl)
    hit = next((i for i, tok in enumerate(tokens[:-1])
                if _program(tok) in BWA_PROGRAMS and tokens[i + 1] == "mem"), None)
    if hit is None:
        return None
    kept = {}
    rest = tokens[hit + 2:]
    j = 0
    while j < len(rest):
        tok = rest[j]
        j += 1
        if not re.match(r"^-[A-Za-z0-9]", tok):
            continue                       # a positional (reference, FASTQ) or a stray fragment
        opt, inline = tok[1], tok[2:]
        if opt in BWA_NUMERIC:
            value = inline
            if not value and j < len(rest):
                value, j = rest[j], j + 1
            kept[opt] = value if NUMERIC_VALUE.match(value or "") else "?"
        elif opt in BWA_FLAGS and not inline:
            kept[opt] = None
        elif opt in BWA_SKIP_VALUE:
            if not inline:
                j += 1                     # skip its value: read group, header text, output path
        # unknown options: arity unknown, value never recorded; the letter is not either
    parts = []
    for opt in sorted(kept, key=lambda o: (o.lower(), o)):
        parts.append(f"-{opt}" if kept[opt] is None else f"-{opt} {kept[opt]}")
    return " ".join(parts) if parts else "default"


def _flag_value(text: str):
    """samtools -F/-f value: decimal, hex, octal, or comma-separated flag names. None if unparsed."""
    text = text.strip()
    try:
        return int(text, 0)
    except ValueError:
        pass
    total = 0
    for name in text.split(","):
        bit = SAMTOOLS_FLAG_NAMES.get(name.strip().upper())
        if bit is None:
            return None
        total |= bit
    return total


# samtools view options that take a value (short, and long without "="): their value is never a
# positional, so it cannot be mistaken for a region.
SAMTOOLS_VIEW_SHORT_VALUE = {"-o", "-U", "-t", "-T", "-L", "-N", "-r", "-R", "-d", "-D", "-q", "-l",
                             "-m", "-x", "-@", "-s", "-e", "-F", "-f", "-G", "-X", "-O"}
SAMTOOLS_VIEW_LONG_VALUE = {
    "--exclude-flags", "--excl-flags", "--require-flags", "--incl-flags", "--rf", "--min-MQ",
    "--min-mq", "--target-file", "--region-file", "--regions-file", "--expr", "--output",
    "--output-fmt", "--output-fmt-option", "--input-fmt-option", "--reference", "--threads",
    "--subsample", "--subsample-seed", "--tag", "--tag-file", "--qname-file", "--library",
    "--read-group", "--read-group-file", "--min-qlen", "--keep-tag", "--remove-tag", "--sanitize",
    "--unoutput", "--fai-reference", "--remove-flags", "--add-flags"}
EXCLUDES_UNMAPPED_EXPR = re.compile(r"!\(?flag\.unmap|\(?flag&4\)?==0|!\(flag&4\)|mapq>0|mapq>=[1-9]")


def _samtools_view_filters(args: list[str]) -> set[str]:
    found, positionals = set(), []
    j = 0
    while j < len(args):
        a = args[j]
        j += 1
        if a.startswith("--"):
            name, eq, value = a.partition("=")
            if not eq and name in SAMTOOLS_VIEW_LONG_VALUE and j < len(args):
                value, j = args[j], j + 1
        elif a.startswith("-") and len(a) > 1:
            name, value = a[:2], a[2:]
            if name in SAMTOOLS_VIEW_SHORT_VALUE and not value and j < len(args):
                value, j = args[j], j + 1
        else:
            positionals.append(a)
            continue
        if name in ("-F", "--exclude-flags", "--excl-flags"):
            bits = _flag_value(value)
            if bits is not None and bits & 0x4:
                found.add("samtools_view_exclude_unmapped")
        elif name in ("-f", "--require-flags"):
            bits = _flag_value(value)
            if bits is not None and bits & 0x2:
                found.add("samtools_view_require_proper_pair")
        elif name in ("-q", "--min-MQ", "--min-mq"):
            if re.fullmatch(r"\d+", value) and int(value) > 0:
                found.add("samtools_view_min_mapq")
        elif name in ("-L", "--target-file", "--region-file", "--regions-file"):
            found.add("samtools_view_region")
        elif name in ("-e", "--expr"):
            if EXCLUDES_UNMAPPED_EXPR.search(value.replace(" ", "")):
                found.add("samtools_view_exclude_unmapped")
    # samtools view IN REGION...: a region list after the input keeps placed reads only.
    if len(positionals) > 1:
        found.add("samtools_view_region")
    return found


def unmapped_filters(cl: str) -> set[str]:
    """FILTER_REASONS found in one @PG CL (read in memory, never returned)."""
    tokens = _split(cl)
    found = set()
    for i, tok in enumerate(tokens[:-1]):
        prog, sub = _program(tok), tokens[i + 1]
        if prog == "samtools" and sub == "view":
            found |= _samtools_view_filters(tokens[i + 2:])
        elif prog == "sambamba" and sub == "view":
            args = tokens[i + 2:]
            for j, a in enumerate(args):
                value = None
                if a in ("-F", "--filter") and j + 1 < len(args):
                    value = args[j + 1]
                elif a.startswith("--filter="):
                    value = a.split("=", 1)[1]
                if value and re.search(r"not\s*\(?\s*unmapped|mapping_quality", value, re.I):
                    found.add("sambamba_view_filter")
    return found


def parse(header_text: str, sample: str) -> dict:
    """The provenance row (PROVENANCE_FIELDS) of one header. Raises ValueError without @SQ."""
    sq, pgs = [], []
    for line in header_text.splitlines():
        if line.startswith("@SQ\t"):
            sq.append(tags(line))
        elif line.startswith("@PG\t"):
            pgs.append(tags(line))
    if not sq:
        raise ValueError("no @SQ lines: an unaligned or truncated header (run unaligned input as FASTQ)")
    names = [t.get("SN", "") for t in sq]
    m5 = [t.get("M5", "").lower() for t in sq]
    chr1 = next((t for t in sq if t.get("SN") in ("chr1", "1")), None)
    chr1_len = int(chr1["LN"]) if chr1 and chr1.get("LN", "").isdigit() else None
    classes = collections.Counter(contig_class(n) for n in names)

    options = [o for o in (bwa_mem_options(p.get("CL", "")) for p in pgs) if o is not None]
    filters = set()
    for p in pgs:
        filters |= unmapped_filters(p.get("CL", ""))
    return {
        "sample": sample,
        "input_ref_fingerprint": md5([f"{t.get('SN', '')}:{t.get('LN', '')}" for t in sq]),
        "input_ref_m5_fingerprint": md5(sorted(set(m5))) if all(m5) else NA,
        "input_ref_contigs": len(sq),
        "input_ref_build": (CHR1_BUILD.get(chr1_len, "unknown") if chr1_len is not None else "no_chr1"),
        "input_ref_alt_contigs": classes.get("alt", 0),
        "input_ref_decoy_contigs": classes.get("decoy", 0),
        "input_ref_hla_contigs": classes.get("hla", 0),
        "input_ref_viral_contigs": classes.get("viral", 0),
        "input_ref_has_ebv": "true" if any(EBV.match(n) for n in names) else "false",
        "input_pg_chain": chain_text(pgs),
        "input_bwa_options": " | ".join(dict.fromkeys(options)) if options else NA,
        "input_unmapped_filter": ",".join(r for r in FILTER_REASONS if r in filters) or "none",
    }


def write_row(path: Path, row: dict) -> None:
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROVENANCE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def summarize(rows: list[dict]) -> dict:
    """Cohort view for RUN_REPORT / the run record: counts only, never a sample ID or a path."""
    def tally(key):
        out = collections.Counter(key(r) for r in rows)
        return dict(sorted(out.items(), key=lambda kv: (-kv[1], kv[0])))

    return {
        "inputs": len(rows),
        "references": tally(lambda r: (
            f"{(r.get('input_ref_fingerprint') or NA)[:12]} build={r.get('input_ref_build') or NA} "
            f"contigs={r.get('input_ref_contigs') or NA} alt={r.get('input_ref_alt_contigs') or NA} "
            f"decoy={r.get('input_ref_decoy_contigs') or NA} hla={r.get('input_ref_hla_contigs') or NA} "
            f"viral={r.get('input_ref_viral_contigs') or NA} ebv={r.get('input_ref_has_ebv') or NA}")),
        "aligners": tally(lambda r: (
            f"{r.get('upstream_aligner') or NA} {r.get('upstream_aligner_version') or NA}"
            + (f" {r['input_bwa_options']}" if r.get("input_bwa_options") not in (None, "", NA) else ""))),
        "pg_chains": tally(lambda r: r.get("input_pg_chain") or NA),
        "unmapped_filter": tally(lambda r: r.get("input_unmapped_filter") or NA),
    }


def main() -> int:
    import argparse
    import sys
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--header", required=True, type=Path, help="SAM header text (samtools view -H)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    try:
        row = parse(a.header.read_text(errors="replace"), a.sample)
    except ValueError as exc:
        sys.exit(f"ERROR: {a.sample}: {exc}")
    write_row(a.output, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
