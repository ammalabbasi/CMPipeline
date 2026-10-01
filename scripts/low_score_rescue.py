#!/usr/bin/env python3
"""Low-score mapped read rescue for BAM/CRAM input (--extract_low_score_mapped auto|on|off).

Why: bwa mem reports a read as unmapped (FLAG 4) only when its best local alignment scores
below -T (default 30). GDC aligned TCGA with `bwa mem -T 0`, so a microbial read with a chance
~19 bp match to human comes out *mapped* (AS < 30, mostly soft-clipped, MAPQ 0-9) and FLAG 4
extraction never sees it. Primary mapped records with AS < 30 are exactly the reads a default
`bwa mem` (-T 30) would have left unmapped. Ported from pksProfiler (feature/e1-extraction,
6138519 / 43eac77), with one difference: `auto` here needs -T < 30 as well as bwa mem, because
with -T >= 30 bwa never reports a primary hit scoring below 30 and the rescue can add nothing.

Two subcommands, both counts/metadata only (no sequence is read or written):

  decide    Reads an alignment header (SAM text from `samtools view -H`), decides per sample
            whether extraction also takes low-score mapped records, prints `on` or `off` on
            stdout, the reason on stderr, and writes one TSV row:
              sample, upstream_aligner, upstream_aligner_version, upstream_bwa_min_score,
              low_as_mode, low_as_mode_reason
            The @PG CL is read in memory (only the program, `mem` and -T are looked at); it is
            never written: command lines carry paths and sample names.

  survival  Counts how many of the rescued read names (the ID list extractReads writes, names
            exactly as `samtools fastq -N` wrote them) are still present after host depletion
            and among the reads KrakenUniq classified.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

# AS:i below this is "low score": bwa mem's default -T. Modules/extract_reads.nf uses the same
# literal; tests/test_low_score_rescue.py fails if the two drift apart.
LOW_AS_THRESHOLD = 30

MODES = ("auto", "on", "off")
# low_as_mode_reason values. build_qc_summary.py and the docs use exactly these.
REASONS = {
    "auto_bwa_mem_T_below_30": "on",      # bwa mem with -T < 30 (GDC TCGA: -T 0)
    "auto_bwa_mem_T_default": "off",      # bwa mem without -T: default 30
    "auto_bwa_mem_T_30_or_more": "off",   # bwa mem with -T >= 30
    "auto_bwa_mem_T_unparsed": "off",     # bwa mem, but -T could not be read
    "auto_bwa_not_mem": "off",            # a bwa @PG that is not `bwa mem` (aln/samse/sampe, no CL)
    "auto_other_aligner": "off",          # STAR, minimap2, bowtie2, ...: AS is not on bwa's scale
    "auto_no_pg": "off",                  # no @PG at all
    "auto_no_aligner_pg": "off",          # only post-processing @PG (samtools sort, Picard, ...)
    "auto_unknown_aligner": "off",        # @PG present, no recognisable aligner
    "requested_on": "on",
    "requested_off": "off",
}
DECISION_FIELDS = ["sample", "upstream_aligner", "upstream_aligner_version",
                   "upstream_bwa_min_score", "low_as_mode", "low_as_mode_reason"]
SURVIVAL_FIELDS = ["sample", "low_as_reads_after_depletion", "low_as_reads_classified"]

BWA_PROGRAMS = ("bwa", "bwa-mem2")
KNOWN_ALIGNERS = ("bwa", "bwa-mem2", "star", "minimap2", "bowtie2", "bowtie", "hisat2",
                  "novoalign", "dragen", "isaac", "tmap", "gsnap", "subread", "bismark",
                  "blasr", "pbmm2", "ngmlr", "segemehl", "smalt", "mosaik")
# Post-processing tools that write @PG but align nothing: never reported as the upstream aligner.
UTILITIES = ("samtools", "picard", "markduplicates", "gatk", "sambamba", "biobambam2", "bamsormadup",
             "bammarkduplicates", "fgbio", "bedtools", "bamtools", "sentieon-util", "umi_tools")
NA = "NA"


def tags(line: str) -> dict:
    out = {}
    for field in line.rstrip("\n").split("\t")[1:]:
        if len(field) > 3 and field[2] == ":":
            out.setdefault(field[:2], field[3:])
    return out


def clean(value: str | None) -> str:
    """A header value safe for a TSV cell and a log line: no tabs/newlines, bounded length."""
    if not value:
        return NA
    return re.sub(r"[^A-Za-z0-9._+:-]", "_", value)[:64] or NA


def _program(token: str) -> str:
    return token.strip("'\"").rsplit("/", 1)[-1]


def _label(pg: dict) -> str:
    """PN, else ID without a numeric / hex de-duplication suffix (bwa.1, bwa-3A2F1C)."""
    return (pg.get("PN") or re.sub(r"([.-]\d+|-[0-9A-F]{6,})$", "", pg.get("ID", ""), flags=re.I)).lower()


def bwa_mem_min_score(cl: str):
    """None if `cl` does not run bwa mem / bwa-mem2 mem; else its -T ('0', 'default', 'unparsed')."""
    tokens = cl.split()
    hit = next((i for i, tok in enumerate(tokens[:-1])
                if _program(tok) in BWA_PROGRAMS and tokens[i + 1] == "mem"), None)
    if hit is None:
        return None
    rest = tokens[hit + 2:]
    for j, tok in enumerate(rest):
        m = re.fullmatch(r"-T(-?\d+)?", tok)
        if m:
            value = m.group(1) if m.group(1) is not None else (rest[j + 1] if j + 1 < len(rest) else "")
            return value if re.fullmatch(r"-?\d+", value) else "unparsed"
    return "default"


def decide(header_text: str, request: str = "auto", sample: str = "") -> dict:
    """The per-sample decision row (see DECISION_FIELDS)."""
    if request not in MODES:
        raise ValueError(f"--extract_low_score_mapped must be auto, on or off (got {request!r})")
    pgs = [tags(line) for line in header_text.splitlines() if line.startswith("@PG\t")]

    bwa_mem = [(pg, bwa_mem_min_score(pg.get("CL", ""))) for pg in pgs]
    bwa_mem = [(pg, t) for pg, t in bwa_mem if t is not None]
    if bwa_mem:
        aligner_pg = bwa_mem[0][0]
        scores = list(dict.fromkeys(t for _, t in bwa_mem))
        numeric = [int(t) for t in scores if re.fullmatch(r"-?\d+", t)]
        if any(t < LOW_AS_THRESHOLD for t in numeric):
            reason = "auto_bwa_mem_T_below_30"
        elif "unparsed" in scores:
            reason = "auto_bwa_mem_T_unparsed"
        elif numeric:
            reason = "auto_bwa_mem_T_30_or_more"
        else:
            reason = "auto_bwa_mem_T_default"
        min_score = ",".join(scores)
        name = aligner_pg.get("PN") or _program(next(
            tok for tok in aligner_pg.get("CL", "").split() if _program(tok) in BWA_PROGRAMS))
    else:
        min_score = NA
        known = [pg for pg in pgs if _label(pg) in KNOWN_ALIGNERS]
        tools = [pg for pg in pgs if _label(pg) not in UTILITIES]
        roots = [pg for pg in tools if "PP" not in pg]
        aligner_pg = known[0] if known else (roots[0] if roots else (tools[0] if tools else None))
        name = (aligner_pg.get("PN") or aligner_pg.get("ID")) if aligner_pg else None
        if not pgs:
            reason = "auto_no_pg"
        elif not tools:
            reason = "auto_no_aligner_pg"
        elif any(_label(pg) in BWA_PROGRAMS for pg in pgs):
            reason = "auto_bwa_not_mem"
        elif known:
            reason = "auto_other_aligner"
        else:
            reason = "auto_unknown_aligner"

    if request != "auto":
        reason = f"requested_{request}"
    return {
        "sample": sample,
        "upstream_aligner": clean(name),
        "upstream_aligner_version": clean(aligner_pg.get("VN")) if aligner_pg else NA,
        "upstream_bwa_min_score": clean(min_score),
        "low_as_mode": REASONS[reason],
        "low_as_mode_reason": reason,
    }


def explain(row: dict) -> str:
    mode = row["low_as_mode"]
    what = (f"also extracting primary mapped records with AS:i < {LOW_AS_THRESHOLD}" if mode == "on"
            else "FLAG 4 records only")
    return (f"low-score mapped extraction {mode} ({row['low_as_mode_reason']}): {what}; "
            f"upstream aligner {row['upstream_aligner']} {row['upstream_aligner_version']}, "
            f"bwa -T {row['upstream_bwa_min_score']}")


def read_names(path: Path):
    """Read names from FASTQ (every 4th line) or FASTA ('>' lines), gzip or plain; first token."""
    with path.open("rb") as probe:
        opener = gzip.open if probe.read(2) == b"\x1f\x8b" else open
    with opener(path, "rt") as handle:
        first = handle.readline()
        if not first:
            return
        if first.startswith(">"):
            yield first[1:].split(maxsplit=1)[0] if first[1:].strip() else ""
            for line in handle:
                if line.startswith(">"):
                    yield line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
            return
        line_no = 0
        line = first
        while line:
            if line_no % 4 == 0:
                if not line.startswith("@"):
                    raise ValueError(f"{path.name}: not FASTQ at record {line_no // 4 + 1}")
                yield line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
            line_no += 1
            line = handle.readline()


def survival(ids: Path, depleted: Path, classified: Path | None, sample: str) -> dict:
    wanted = set()
    with gzip.open(ids, "rt") as handle:
        for line in handle:
            name = line.rstrip("\n")
            if name:
                wanted.add(name)
    after = sum(1 for name in read_names(depleted) if name in wanted)
    if classified is None:
        hit = NA
    else:
        hit = sum(1 for name in read_names(classified) if name in wanted)
    return {"sample": sample, "low_as_reads_after_depletion": after, "low_as_reads_classified": hit}


def write_row(path: Path, fields: list[str], row: dict) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("decide")
    d.add_argument("--header", required=True, type=Path, help="SAM header text (samtools view -H)")
    d.add_argument("--sample", required=True)
    d.add_argument("--mode", default="auto", help="--extract_low_score_mapped: auto | on | off")
    d.add_argument("--output", required=True, type=Path)
    s = sub.add_parser("survival")
    s.add_argument("--ids", required=True, type=Path, help="gzipped rescued read names, one per line")
    s.add_argument("--depleted", required=True, type=Path, help="host-depleted FASTQ(.gz)")
    s.add_argument("--classified", type=Path, default=None, help="KrakenUniq classified FASTA")
    s.add_argument("--sample", required=True)
    s.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()

    if a.cmd == "decide":
        try:
            row = decide(a.header.read_text(errors="replace"), a.mode, a.sample)
        except ValueError as exc:
            sys.exit(f"ERROR: {a.sample}: {exc}")
        write_row(a.output, DECISION_FIELDS, row)
        print(f"extractReads: {a.sample}: {explain(row)}", file=sys.stderr)
        print(row["low_as_mode"])
        return 0

    row = survival(a.ids, a.depleted, a.classified, a.sample)
    write_row(a.output, SURVIVAL_FIELDS, row)
    print(f"lowAsSurvival: {a.sample}: {row['low_as_reads_after_depletion']} rescued reads after "
          f"host depletion, {row['low_as_reads_classified']} classified", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
