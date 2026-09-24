"""Shared helpers for alignment validation (validation_*.py).

Reference contigs are named  <label>|<accession>|<contig>, where <label> is the Kraken species
name with spaces replaced by underscores (so it matches the Bracken tables exactly), or
"Homo_sapiens" for the competitive human genome.
"""
import re

HUMAN_LABEL = "Homo_sapiens"
BRACKEN_SUFFIX = re.compile(r"\.bracken\.[GS](\.mpa)?\.krakenreport\.txt$")


def label_of(species_name):
    """Reference label for a species name: 'Escherichia coli' -> 'Escherichia_coli'.

    Idempotent on already-underscored names. Bracken tables made by kreport2mpa.py replace
    spaces with underscores ('s__Escherichia_coli') while taxDB and RefSeq use spaces, so every
    comparison between the two MUST go through this function (audit B01, 2026-09-23)."""
    return re.sub(r"[\s_]+", "_", species_name.strip())


def spaced(species_name):
    """'Escherichia_coli' -> 'Escherichia coli' (the taxDB / RefSeq spelling)."""
    return re.sub(r"[\s_]+", " ", species_name.strip())


def species_from_clade(clade_name):
    """'d__Bacteria|...|g__Escherichia|s__Escherichia coli' -> 'Escherichia coli' (or None)."""
    m = re.search(r"\|s__([^|]+)$", clade_name)
    return m.group(1) if m else None


def genus_clade(clade_name):
    """Species clade -> its genus clade: drop the trailing |s__... level."""
    return re.sub(r"\|s__[^|]+$", "", clade_name)


def sample_name(column):
    """Bracken/consensus column name -> sample ID (strips a Bracken file suffix if present)."""
    return BRACKEN_SUFFIX.sub("", column)


def read_species_table(path):
    """Merged species table -> (sample IDs, {clade_name: {sample: count}}).

    Accepts the consensus output (first column clade_name) or a raw combine_mpa table
    (first column #Classification). Rows without an s__ level are ignored.
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = [sample_name(c) for c in header[1:]]
        rows = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or species_from_clade(f[0]) is None:
                continue
            rows[f[0]] = {s: float(v or 0) for s, v in zip(samples, f[1:])}
    return samples, rows


def read_groups(path):
    """Groups TSV (group<TAB>member_species_name) -> {member_label: group_label}."""
    groups = {}
    if not path:
        return groups
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            group, member = line.rstrip("\n").split("\t")[:2]
            if group == "group":          # header
                continue
            groups[label_of(member)] = label_of(group)
    return groups
