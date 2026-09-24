#!/usr/bin/env python3
"""Cohort QC summary: one row per sample, from artefacts the pipeline already produces.

Answers the questions a run leaves open:

  - how many reads survived each stage, and what fraction of the input is left
  - did this sample clear --consensus_min_reads, and by how much (so "why has this
    sample no MetaPhlAn profile?" has an answer)
  - how many genera and species each profiler called, and how many they agree on
  - what role decontamination gave the sample

Inputs are globbed from a directory of staged files, so a stage that did not run simply
contributes nothing and its columns read NA. Standard library only.

Counts and rates only -- no sequence, and no clinical value, is read or written.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

NA = "NA"

COLUMNS = [
    "sample",
    # library denominator, BAM/CRAM only (audit C07 / pksProfiler F09). For BAM/CRAM input,
    # reads_raw below is the UNMAPPED subset the pipeline extracted, not the library; these
    # columns carry the library itself. NA for FASTQ input, where reads_raw is the library.
    "library_primary_records", "extracted_unmapped_records", "pct_library_unmapped",
    # input and filtering (fastp)
    "reads_raw", "bases_raw", "q30_rate_raw",
    "reads_after_fastp", "q30_rate_filtered", "duplication_rate",
    "reads_low_quality", "reads_too_short",
    "adapter_scan_applied", "adapter_trimmed_reads", "pct_pass_fastp",
    # host depletion
    "reads_after_hg38", "reads_after_t2t_phix", "reads_host_depleted",
    "pct_removed_hg38", "pct_removed_t2t_phix", "pct_removed_pangenome",
    "pct_retained_overall",      # host-depleted / reads entering fastp (the extracted subset for BAM/CRAM)
    "pct_retained_of_library",   # host-depleted / primary records in the BAM/CRAM; NA for FASTQ
    # classification
    "bracken_microbial_reads", "cleared_consensus_gate",
    "bracken_genera", "bracken_species",
    "metaphlan_genera", "metaphlan_species",
    "shared_genera", "shared_species",
    # decontamination
    "decontam_sample_type", "decontam_role", "decontam_final_status",
    # what happened to the sample, in one word (audit C13/C02/C10): no_outputs (in the sample sheet
    # but nothing reached QC: a failed or --sample_failure_strategy ignore'd task),
    # dropped_zero_after_fastp, no_reads_after_depletion, below_bracken_threshold (Bracken wrote a
    # header-only no-call), profiled; NA when the stages that decide it did not run
    "bracken_status", "sample_status",
]


def pct(numerator, denominator):
    # A zero count over a real denominator is 0.00, not NA (audit C13): NA only when a value is
    # genuinely unknown or the denominator is zero.
    if numerator in (None, "", NA) or denominator in (None, "", NA):
        return NA
    try:
        return f"{100.0 * float(numerator) / float(denominator):.2f}"
    except (TypeError, ValueError, ZeroDivisionError):
        return NA


def removed(before, after):
    if before in (None, NA) or after in (None, NA):
        return NA
    try:
        before, after = float(before), float(after)
    except (TypeError, ValueError):
        return NA
    if before <= 0:
        return NA
    return f"{100.0 * (before - after) / before:.2f}"


# ------------------------------------------------------------------------ parsers


def parse_fastp(path: Path) -> dict:
    data = json.loads(path.read_text())
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    filtered = data.get("filtering_result", {})
    adapter = data.get("adapter_cutting")
    row = {
        "reads_raw": before.get("total_reads"),
        "bases_raw": before.get("total_bases"),
        "q30_rate_raw": before.get("q30_rate"),
        "reads_after_fastp": after.get("total_reads"),
        "q30_rate_filtered": after.get("q30_rate"),
        "duplication_rate": data.get("duplication", {}).get("rate"),
        "reads_low_quality": filtered.get("low_quality_reads"),
        "reads_too_short": filtered.get("too_short_reads"),
        # `adapter_cutting` appears only when fastp was given --adapter_fasta, which is
        # exactly the per-sample decision --adapter_trim auto makes.
        "adapter_scan_applied": "true" if adapter else "false",
        "adapter_trimmed_reads": (adapter or {}).get("adapter_trimmed_reads", NA),
    }
    row["pct_pass_fastp"] = pct(row["reads_after_fastp"], row["reads_raw"])
    return row


def parse_extract_counts(path: Path) -> dict:
    with path.open() as handle:
        record = next(csv.DictReader(handle, delimiter="\t"), None) or {}
    return {
        "library_primary_records": record.get("library_primary_records") or NA,
        "extracted_unmapped_records": record.get("extracted_unmapped_records") or NA,
    }


def parse_depletion(path: Path) -> dict:
    stages = {}
    with path.open() as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            stages[record["stage"]] = record["reads"]
    return {
        "reads_after_hg38": stages.get("after_hg38", NA),
        "reads_after_t2t_phix": stages.get("after_t2t_phix", NA),
        "reads_host_depleted": stages.get("host_depleted", NA),
        "_depletion_input": stages.get("input"),
    }


def parse_krakenuniq_microbial_reads(path: Path):
    """Reads assigned to any taxon -- the identical definition Bracken.nf gates on.

    Keeping one definition matters: a QC table that disagreed with the gate would be
    worse than no table.

    A KrakenUniq report has 9 tab-separated columns:
        % reads taxReads kmers dup cov taxID rank taxName
    Both the root and unclassified rows carry rank "no rank", so they are identified by
    taxName, not by rank. Root and unclassified are siblings: root's own clade count is
    the classified total, and subtracting unclassified from it is wrong.
    """
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[8].strip() == "root":
                try:
                    return int(float(fields[1].strip()))
                except ValueError:
                    return None
    return None


def parse_bracken_taxa(path: Path) -> set[str]:
    """Taxon names with a non-zero estimate."""
    names = set()
    with path.open() as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            name = (record.get("name") or "").strip()
            estimate = record.get("new_est_reads") or record.get("kraken_assigned_reads") or "0"
            try:
                if name and float(estimate) > 0:
                    names.add(name)
            except ValueError:
                continue
    return names


def parse_metaphlan_taxa(path: Path) -> tuple[set[str], set[str]]:
    """Genus and species names from a MetaPhlAn profile.

    Names are normalised to Bracken's spelling -- `g__Escherichia` -> `Escherichia`,
    `s__Escherichia_coli` -> `Escherichia coli` -- so the two can be intersected.
    """
    genera, species = set(), set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            clade = line.split("\t")[0]
            if "|s__" in clade or clade.startswith("s__"):
                species.add(clade.rsplit("s__", 1)[1].replace("_", " ").strip())
            elif "|g__" in clade or clade.startswith("g__"):
                genera.add(clade.rsplit("g__", 1)[1].replace("_", " ").strip())
    return genera, species


def parse_sample_validation(path: Path) -> dict[str, dict]:
    roles = {}
    with path.open() as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            sample = (record.get("sample_id") or "").strip()
            if sample:
                roles[sample] = {
                    "decontam_sample_type": record.get("sample_type") or NA,
                    "decontam_role": record.get("role") or NA,
                    "decontam_final_status": record.get("final_status") or NA,
                }
    return roles


def sample_status(row: dict) -> str:
    def is_zero(key):
        return str(row.get(key)) in ("0", "0.0")
    if is_zero("reads_after_fastp"):
        return "dropped_zero_after_fastp"
    if is_zero("reads_host_depleted"):
        return "no_reads_after_depletion"
    if row.get("bracken_status") == "no_call":
        return "below_bracken_threshold"
    if row.get("bracken_status") == "called":
        return "profiled"
    return NA


# --------------------------------------------------------------------------- build


def sample_of(path: Path, suffix: str) -> str:
    return path.name[: -len(suffix)] if path.name.endswith(suffix) else path.stem


def build(inputs: Path, min_reads: int) -> list[dict]:
    rows: dict[str, dict] = {}

    def row_for(sample: str) -> dict:
        return rows.setdefault(sample, {c: NA for c in COLUMNS} | {"sample": sample})

    for path in sorted(inputs.rglob("*.fastp.json")):
        row_for(sample_of(path, ".fastp.json")).update(parse_fastp(path))

    for path in sorted(inputs.rglob("*.extract_counts.tsv")):
        row_for(sample_of(path, ".extract_counts.tsv")).update(parse_extract_counts(path))

    for path in sorted(inputs.rglob("*.depletion_counts.tsv")):
        row_for(sample_of(path, ".depletion_counts.tsv")).update(parse_depletion(path))

    for path in sorted(inputs.rglob("*.krakenuniq.report.txt")):
        sample = sample_of(path, ".krakenuniq.report.txt")
        reads = parse_krakenuniq_microbial_reads(path)
        row = row_for(sample)
        row["bracken_microbial_reads"] = NA if reads is None else reads
        row["cleared_consensus_gate"] = (
            NA if reads is None else ("true" if reads >= min_reads else "false")
        )

    bracken_genera: dict[str, set[str]] = {}
    bracken_species: dict[str, set[str]] = {}
    for path in sorted(inputs.rglob("*.bracken.G.report.txt")):
        sample = sample_of(path, ".bracken.G.report.txt")
        bracken_genera[sample] = parse_bracken_taxa(path)
        row = row_for(sample)
        row["bracken_genera"] = len(bracken_genera[sample])
        # Bracken.nf writes a header-only report when a sample is below its threshold (C02): that is
        # a no-call, not "ran and found 0 taxa".
        with path.open() as handle:
            data_rows = sum(1 for line in handle if line.strip()) - 1
        row["bracken_status"] = "called" if data_rows > 0 else "no_call"
    for path in sorted(inputs.rglob("*.bracken.S.report.txt")):
        sample = sample_of(path, ".bracken.S.report.txt")
        bracken_species[sample] = parse_bracken_taxa(path)
        row_for(sample)["bracken_species"] = len(bracken_species[sample])

    for path in sorted(inputs.rglob("*.profiled_metagenome.txt")):
        sample = sample_of(path, ".profiled_metagenome.txt")
        genera, species = parse_metaphlan_taxa(path)
        row = row_for(sample)
        row["metaphlan_genera"] = len(genera)
        row["metaphlan_species"] = len(species)
        # Agreement between the two profilers, per sample. The consensus step computes a
        # cohort-level whitelist; this is the per-sample view, which the cohort figure
        # cannot show.
        if sample in bracken_genera:
            row["shared_genera"] = len(genera & bracken_genera[sample])
        if sample in bracken_species:
            row["shared_species"] = len(species & bracken_species[sample])

    for path in sorted(inputs.rglob("*sample_validation.tsv")):
        for sample, roles in parse_sample_validation(path).items():
            if sample in rows:
                rows[sample].update(roles)

    # Samples the sheet asked for that produced nothing at all (audit C10).
    for path in sorted(inputs.rglob("expected_samples.txt")):
        for sample in (s.strip() for s in path.read_text().splitlines()):
            if sample and sample not in rows:
                row_for(sample)["sample_status"] = "no_outputs"

    # derived retention, once every count is in
    for row in rows.values():
        # A sample fastp emptied never reaches depletion: its depletion counts are 0, not unknown.
        if str(row.get("reads_after_fastp")) == "0":
            for column in ("reads_after_hg38", "reads_after_t2t_phix", "reads_host_depleted"):
                if row.get(column) in (None, NA):
                    row[column] = 0
        if row.get("sample_status") != "no_outputs":
            row["sample_status"] = sample_status(row)
        depletion_input = row.pop("_depletion_input", None)
        first = depletion_input or row.get("reads_after_fastp")
        row["pct_removed_hg38"] = removed(first, row.get("reads_after_hg38"))
        row["pct_removed_t2t_phix"] = removed(row.get("reads_after_hg38"),
                                              row.get("reads_after_t2t_phix"))
        row["pct_removed_pangenome"] = removed(row.get("reads_after_t2t_phix"),
                                               row.get("reads_host_depleted"))
        row["pct_retained_overall"] = pct(row.get("reads_host_depleted"), row.get("reads_raw"))
        row["pct_library_unmapped"] = pct(row.get("extracted_unmapped_records"),
                                          row.get("library_primary_records"))
        row["pct_retained_of_library"] = pct(row.get("reads_host_depleted"),
                                             row.get("library_primary_records"))

    return [rows[s] for s in sorted(rows)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path,
                        help="directory of staged per-sample QC artefacts")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-reads", type=int, default=100000,
                        help="--consensus_min_reads, for cleared_consensus_gate")
    args = parser.parse_args()

    rows = build(args.inputs, args.min_reads)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"QC summary: {len(rows)} sample(s) -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
