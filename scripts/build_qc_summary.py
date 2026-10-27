#!/usr/bin/env python3
"""Cohort QC summary: one row per sample, from artefacts the pipeline already produces.

Answers the questions a run leaves open:

  - how many reads survived each stage, and what fraction of the input is left
  - did this sample clear --consensus_min_reads, and by how much (so "why has this
    sample no MetaPhlAn profile?" has an answer)
  - how many genera and species each profiler called, and how many they agree on
  - what role decontamination gave the sample
  - whether batch correction included it, per corrected rank

Inputs are globbed from a directory of staged files, so a stage that did not run simply
contributes nothing and its columns read NA. A run that starts at decontamination or batch
correction gets its rows from those stages. Standard library only.

Counts and rates only -- no sequence, and no clinical value, is read or written.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
from pathlib import Path

NA = "NA"

COLUMNS = [
    "sample",
    # library denominator, BAM/CRAM only (audit C07 / pksProfiler F09). For BAM/CRAM input,
    # reads_raw below is the UNMAPPED subset the pipeline extracted, not the library; these
    # columns carry the library itself. NA for FASTQ input, where reads_raw is the library.
    "library_primary_records", "extracted_unmapped_records", "pct_library_unmapped",
    # low-score mapped rescue, BAM/CRAM only (--extract_low_score_mapped; Ludmil 1 Oct, L8). The
    # upstream aligner, version and bwa -T come from @PG; low_as_mode is on/off and
    # low_as_mode_reason why (scripts/low_score_rescue.py REASONS). low_as_reads_extracted are
    # primary MAPPED records with AS < 30 added to the FLAG 4 reads (so reads_raw = FLAG 4 +
    # rescued); then how many of them survive host depletion and how many KrakenUniq classified.
    # low_as_share_of_bracken = low_as_reads_classified / kraken_classified_reads, in percent: an
    # estimate of the rescued share of the table Bracken re-estimates from (Bracken itself does not
    # attribute reads). NA for FASTQ input; the read counts are 0 when the rescue was off
    # (low_as_no_as_records is NA then: records without AS are counted only when it is on).
    "upstream_aligner", "upstream_aligner_version", "upstream_bwa_min_score",
    "low_as_mode", "low_as_mode_reason", "low_as_reads_extracted", "low_as_no_as_records",
    "low_as_reads_after_depletion", "low_as_reads_classified", "low_as_share_of_bracken",
    # input provenance, BAM/CRAM only (audit X13; scripts/header_provenance.py, from the header the
    # rescue decision reads): the @SQ fingerprint (md5 of ordered SN:LN) and M5-set fingerprint,
    # contig count, build from chr1's length, contig classes (ALT / hs38d1 decoy / HLA / viral incl.
    # chrEBV), the @PG chain as program:version, bwa mem options reduced to letters and numbers, and
    # any upstream step that removed unmapped reads (none | FILTER_REASONS). NA for FASTQ input.
    "input_ref_fingerprint", "input_ref_m5_fingerprint", "input_ref_contigs", "input_ref_build",
    "input_ref_alt_contigs", "input_ref_decoy_contigs", "input_ref_hla_contigs",
    "input_ref_viral_contigs", "input_ref_has_ebv", "input_pg_chain", "input_bwa_options",
    "input_unmapped_filter",
    # input and filtering (fastp)
    "reads_raw", "bases_raw", "q30_rate_raw",
    "reads_after_fastp", "q30_rate_filtered", "duplication_rate",
    "reads_low_quality", "reads_too_short",
    "adapter_scan_applied", "adapter_list", "adapter_trimmed_reads", "pct_pass_fastp",
    # host depletion
    "reads_after_hg38", "reads_after_t2t_phix", "reads_host_depleted",
    "pct_removed_hg38", "pct_removed_t2t_phix", "pct_removed_pangenome",
    "pct_retained_overall",      # host-depleted / reads entering fastp (the extracted subset for BAM/CRAM)
    "pct_retained_of_library",   # host-depleted / primary records in the BAM/CRAM; NA for FASTQ
    # classification
    # KrakenUniq clade counts: classified = root (what --consensus_min_reads gates on), which
    # includes human (9606) and synthetic construct (32630, vectors/adapters); microbial excludes both.
    "kraken_classified_reads", "kraken_human_reads", "kraken_synthetic_reads",
    "bracken_microbial_reads", "cleared_consensus_gate",
    # Bracken species estimates for the same two taxa, and their share of all species reads
    "bracken_human_reads", "bracken_synthetic_reads", "pct_bracken_human", "pct_bracken_synthetic",
    "bracken_genera", "bracken_species",
    "metaphlan_genera", "metaphlan_species",
    "shared_genera", "shared_species",
    # decontamination
    "decontam_sample_type", "decontam_role", "decontam_final_status",
    # batch correction, per corrected rank: included, or excluded: <reason>
    "batch_species_status", "batch_genus_status",
    # what happened to the sample, in one word (audit C13/C02/C10): no_outputs (in the sample sheet
    # but nothing reached QC: a failed or --sample_failure_strategy ignore'd task),
    # dropped_zero_after_fastp, no_reads_after_depletion, below_bracken_threshold (Bracken wrote a
    # header-only no-call), profiled; NA when the stages that decide it did not run
    "bracken_status", "sample_status",
    # Gated per-sample stages (audit X05): ran | below_gate | failed | not_requested, and for HUMAnN
    # also metaphlan_failed. `failed` means the sample cleared --consensus_min_reads (and, for HUMAnN,
    # MetaPhlAn ran) but the stage left no output: a failed task under --sample_failure_strategy
    # ignore, never "no microbes". NA when classification itself did not finish for the sample.
    "metaphlan_status", "humann_status",
    # QC floor (pksProfiler COHORT_CHECKS A3; port review #34): pass, or why the sample cannot be
    # assessed -- unmapped_filtered_upstream (input_unmapped_filter is not none), no_unmapped_records
    # (zero FLAG 4 records in the BAM/CRAM), below_min_reads_after_depletion (fewer host-depleted
    # reads than --min_reads_after_depletion, when > 0). NA when not evaluated. A not-assessable
    # sample is never a microbial negative: terminal_state not_assessable_<reason>.
    "qc_floor",
    # How the sample ENDED, exactly one of (audit X05/X06; pksProfiler e126439):
    #   completed_profiled | completed_no_call | completed   every per-sample stage it needed finished
    #   failed_at_<stage>        the first stage it needed that left no output (a failed task)
    #   excluded_<reason>        dropped by the pipeline on purpose: zero_reads_after_fastp, or
    #                            decontam_<reason> (decontamination's sample_validation.tsv)
    #   not_assessable_<reason>  the QC floor (qc_floor): the input cannot say whether microbes are
    #                            present -- never a negative. Takes precedence over every state but
    #                            failed_at_extraction / earlier stage failures it depends on
    # Stages, in order: extraction (BAM/CRAM only), read_qc, host_depletion, classification,
    # metaphlan (cleared samples only), humann (when requested, after MetaPhlAn ran),
    # low_as_survival (rescue-on samples only). See terminal_state().
    "terminal_state",
]

# Per-sample stages in pipeline order (--stages takes a subset; main.nf passes the ones this run
# requested). A sample that needed a stage and has no output from it failed there; stages after a
# failure never run, so the FIRST gap is the failure.
STAGE_ORDER = ("extraction", "read_qc", "host_depletion", "classification", "metaphlan", "humann",
               "low_as_survival")
ALIGNMENT_INPUTS = {"bam", "cram"}


def _known(value) -> bool:
    return value not in (None, "", NA)


def stage_needed(stage: str, row: dict, input_type: str | None) -> bool:
    """Whether this sample had to pass through `stage`, given what earlier stages decided."""
    if stage == "extraction":
        # Only BAM/CRAM inputs are extracted; without the input type (a hand run) it is not checked.
        return input_type in ALIGNMENT_INPUTS
    if stage == "metaphlan":
        return row.get("cleared_consensus_gate") == "true"
    if stage == "humann":
        return row.get("cleared_consensus_gate") == "true" and _known(row.get("metaphlan_genera"))
    if stage == "low_as_survival":
        return row.get("low_as_mode") == "on"
    return True


def stage_done(stage: str, row: dict, humann_done: set) -> bool:
    """Whether the stage's per-sample output reached the QC summary."""
    return {
        "extraction": lambda: _known(row.get("library_primary_records")),
        "read_qc": lambda: _known(row.get("reads_after_fastp")),
        "host_depletion": lambda: _known(row.get("reads_host_depleted")),
        "classification": lambda: _known(row.get("bracken_status")) and _known(row.get("bracken_microbial_reads")),
        "metaphlan": lambda: _known(row.get("metaphlan_genera")),
        "humann": lambda: row["sample"] in humann_done,
        "low_as_survival": lambda: _known(row.get("low_as_reads_after_depletion")),
    }[stage]()


def gated_status(stage: str, row: dict, stages: list[str], humann_done: set) -> str:
    """metaphlan_status / humann_status: ran | below_gate | failed | metaphlan_failed | not_requested | NA."""
    if stage not in stages:
        return "not_requested"
    gate = row.get("cleared_consensus_gate")
    if gate == "false":
        return "below_gate"
    if gate != "true":
        return NA
    if stage == "humann" and not _known(row.get("metaphlan_genera")):
        return "metaphlan_failed" if "metaphlan" in stages else NA
    return "ran" if stage_done(stage, row, humann_done) else "failed"


QC_FLOOR_INPUT = ("unmapped_filtered_upstream", "no_unmapped_records")
QC_FLOOR_DEPTH = "below_min_reads_after_depletion"


def _int(value):
    try:
        return None if value in (None, "", NA) else int(float(value))
    except (TypeError, ValueError):
        return None


def qc_floor(row: dict, input_type: str | None, min_after_depletion: int = 0) -> str:
    """The QC floor (see COLUMNS): pass | unmapped_filtered_upstream | no_unmapped_records |
    below_min_reads_after_depletion | NA (nothing to judge it on)."""
    alignment = input_type in ALIGNMENT_INPUTS or (input_type is None and _known(row.get("library_primary_records")))
    judged = False
    if alignment:
        upstream = row.get("input_unmapped_filter")
        if _known(upstream):
            judged = True
            if upstream != "none":
                return "unmapped_filtered_upstream"
        unmapped = _int(row.get("extracted_unmapped_records"))
        if unmapped is not None:
            judged = True
            if unmapped == 0:
                return "no_unmapped_records"
    # A sample fastp emptied is excluded_zero_reads_after_fastp (unchanged); the depth floor is for
    # samples that reached host depletion with reads.
    if min_after_depletion > 0 and str(row.get("reads_after_fastp")) not in ("0", "0.0"):
        depleted = _int(row.get("reads_host_depleted"))
        if depleted is not None:
            judged = True
            if depleted < min_after_depletion:
                return QC_FLOOR_DEPTH
    return "pass" if judged else NA


def terminal_state(row: dict, stages: list[str], input_type: str | None, humann_done: set,
                   per_sample: bool) -> str:
    """Exactly one terminal state per sample (see COLUMNS)."""
    if per_sample:
        floor = row.get("qc_floor")
        # An input-level floor is known as soon as extraction finished; everything after it was
        # skipped on purpose (--qc_floor_action exclude), so it is not a failure.
        if floor in QC_FLOOR_INPUT:
            return f"not_assessable_{floor}"
        if str(row.get("reads_after_fastp")) == "0":
            return "excluded_zero_reads_after_fastp"
        if floor == QC_FLOOR_DEPTH:
            return f"not_assessable_{floor}"
        for stage in stages:
            if stage_needed(stage, row, input_type) and not stage_done(stage, row, humann_done):
                return f"failed_at_{stage}"
    decontam = row.get("decontam_final_status")
    if _known(decontam) and decontam != "included":
        reason = re.sub(r"[^a-z0-9]+", "_", (row.get("_decontam_reason") or decontam).lower()).strip("_")
        return f"excluded_decontam_{reason or 'unspecified'}"
    if row.get("bracken_status") == "called":
        return "completed_profiled"
    if row.get("bracken_status") == "no_call":
        return "completed_no_call"
    return "completed"


def parse_stages(text: str | None) -> list[str]:
    stages = [s.strip() for s in (text or "").split(",") if s.strip()]
    unknown = [s for s in stages if s not in STAGE_ORDER]
    if unknown:
        raise ValueError(f"Unknown stage(s) {', '.join(unknown)}; expected: {', '.join(STAGE_ORDER)}")
    return [s for s in STAGE_ORDER if s in stages]


def infer_stages(rows: dict) -> list[str]:
    """For hand runs without --stages: every stage at least one sample reported output for. A stage
    that every sample failed cannot be inferred, which is why main.nf passes --stages."""
    seen = lambda col: any(_known(r.get(col)) for r in rows.values())
    stages = []
    for stage, col in (("extraction", "library_primary_records"), ("read_qc", "reads_after_fastp"),
                       ("host_depletion", "reads_host_depleted"), ("classification", "bracken_status"),
                       ("metaphlan", "metaphlan_genera"), ("low_as_survival", "low_as_reads_after_depletion")):
        if seen(col):
            stages.append(stage)
    return stages


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


def fastp_adapter_list(command: str) -> str:
    """Basename of the --adapter_fasta fastp ran with, `none`, or NA without a command line."""
    if not command:
        return NA
    words = shlex.split(command)
    for i, word in enumerate(words):
        if word.startswith("--adapter_fasta="):
            return Path(word.split("=", 1)[1]).name
        if word == "--adapter_fasta" and i + 1 < len(words):
            return Path(words[i + 1]).name
    return "none"


def parse_fastp(path: Path) -> dict:
    data = json.loads(path.read_text())
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    filtered = data.get("filtering_result", {})
    adapter = data.get("adapter_cutting")
    adapter_list = fastp_adapter_list(data.get("command", ""))
    row = {
        "reads_raw": before.get("total_reads"),
        "bases_raw": before.get("total_bases"),
        "q30_rate_raw": before.get("q30_rate"),
        "reads_after_fastp": after.get("total_reads"),
        "q30_rate_filtered": after.get("q30_rate"),
        "duplication_rate": data.get("duplication", {}).get("rate"),
        "reads_low_quality": filtered.get("low_quality_reads"),
        "reads_too_short": filtered.get("too_short_reads"),
        # Which adapter FASTA fastp was given, read from its own command line. fastp writes
        # `adapter_cutting` whenever it trimmed adapters, including by auto-detection alone
        # (seen with 0.24 and 1.3.7), so that block cannot tell whether a FASTA was scanned.
        # adapter_trimmed_reads counts both.
        "adapter_scan_applied": NA if adapter_list == NA else ("false" if adapter_list == "none" else "true"),
        "adapter_list": adapter_list,
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
        # Absent from extract_counts.tsv written before the low-score rescue existed: NA.
        "low_as_reads_extracted": record.get("low_as_reads_extracted") or NA,
        "low_as_no_as_records": record.get("low_as_no_as_records") or NA,
        # What extractReads decided (QC floor); cross-checked against the recomputed qc_floor.
        "_input_qc_floor": record.get("input_qc_floor") or None,
    }


LOW_AS_MODES = {"on", "off"}
UPSTREAM_FIELDS = ("upstream_aligner", "upstream_aligner_version", "upstream_bwa_min_score",
                   "low_as_mode", "low_as_mode_reason")


def parse_upstream_aligner(path: Path) -> dict:
    """<sample>.upstream_aligner.tsv from extractReads (scripts/low_score_rescue.py decide)."""
    with path.open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != 1:
        raise ValueError(f"{path.name}: expected one row, got {len(records)}")
    record = records[0]
    if record.get("low_as_mode") not in LOW_AS_MODES:
        raise ValueError(f"{path.name}: low_as_mode must be on or off, got {record.get('low_as_mode')!r}")
    return {field: record.get(field) or NA for field in UPSTREAM_FIELDS}


PROVENANCE_FIELDS = ("input_ref_fingerprint", "input_ref_m5_fingerprint", "input_ref_contigs",
                     "input_ref_build", "input_ref_alt_contigs", "input_ref_decoy_contigs",
                     "input_ref_hla_contigs", "input_ref_viral_contigs", "input_ref_has_ebv",
                     "input_pg_chain", "input_bwa_options", "input_unmapped_filter")


def parse_input_provenance(path: Path) -> dict:
    """<sample>.input_provenance.tsv from extractReads (scripts/header_provenance.py)."""
    with path.open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != 1:
        raise ValueError(f"{path.name}: expected one row, got {len(records)}")
    return {field: records[0].get(field) or NA for field in PROVENANCE_FIELDS}


def parse_low_as_survival(path: Path) -> dict:
    with path.open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != 1:
        raise ValueError(f"{path.name}: expected one row, got {len(records)}")
    record = records[0]
    return {field: record.get(field) or NA
            for field in ("low_as_reads_after_depletion", "low_as_reads_classified")}


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


HUMAN_TAXID = "9606"
SYNTHETIC_TAXID = "32630"   # "synthetic construct": UniVec/EmVec vectors, adapters, primers


def parse_krakenuniq_clades(path: Path) -> dict:
    """Clade read counts for root, human and synthetic construct, keyed by column name.

    Clade count = column 2 (taxReads) of the 9-column report; taxID is column 7. A taxon
    the report does not list had no reads (0). Root is required: without it every value is NA.
    """
    counts = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            taxid, name = fields[6].strip(), fields[8].strip()
            key = ("root" if name == "root" else "human" if taxid == HUMAN_TAXID
                   else "synthetic" if taxid == SYNTHETIC_TAXID else None)
            if key and key not in counts:
                try:
                    counts[key] = int(float(fields[1].strip()))
                except ValueError:
                    counts[key] = None
    root = counts.get("root")
    if root is None:
        return {c: NA for c in ("kraken_classified_reads", "kraken_human_reads",
                                "kraken_synthetic_reads", "bracken_microbial_reads")}
    human, synthetic = counts.get("human") or 0, counts.get("synthetic") or 0
    return {"kraken_classified_reads": root, "kraken_human_reads": human,
            "kraken_synthetic_reads": synthetic,
            "bracken_microbial_reads": max(root - human - synthetic, 0)}


def parse_bracken_nonmicrobial(path: Path) -> dict:
    """Bracken species estimates for human and synthetic construct, and their % of all species reads."""
    total, found = 0.0, {HUMAN_TAXID: 0.0, SYNTHETIC_TAXID: 0.0}
    with path.open() as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            try:
                reads = float(record.get("new_est_reads") or 0)
            except ValueError:
                continue
            total += reads
            taxid = (record.get("taxonomy_id") or "").strip()
            if taxid in found:
                found[taxid] += reads
    human, synthetic = int(found[HUMAN_TAXID]), int(found[SYNTHETIC_TAXID])
    return {"bracken_human_reads": human, "bracken_synthetic_reads": synthetic,
            "pct_bracken_human": pct(human, int(total)) if total else NA,
            "pct_bracken_synthetic": pct(synthetic, int(total)) if total else NA}


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
        reader = csv.DictReader(handle, delimiter="\t")
        if "in_otu" not in (reader.fieldnames or []) or "final_status" not in (reader.fieldnames or []):
            raise ValueError(f"{path.name}: sample_validation needs in_otu and final_status columns")
        for record in reader:
            sample = (record.get("sample_id") or "").strip()
            if not sample:
                raise ValueError(f"{path.name}: row without a sample_id")
            in_otu = record["in_otu"]
            if in_otu not in ("TRUE", "FALSE"):
                raise ValueError(f"{path.name}: in_otu must be TRUE or FALSE for {sample}, got {in_otu!r}")
            if sample in roles:
                raise ValueError(f"{path.name}: duplicate sample {sample}")
            roles[sample] = {
                "decontam_sample_type": record.get("sample_type") or NA,
                "decontam_role": record.get("role") or NA,
                "decontam_final_status": record.get("final_status") or NA,
                # Why an excluded sample was excluded (metadata_only, insufficient_library, ...);
                # used for terminal_state only, not a column.
                "_decontam_reason": record.get("reason") or "",
                # FALSE for a metadata-only row: the sample is not in the analysed table.
                "_in_table": in_otu == "TRUE",
            }
    return roles


def parse_batch_sample_status(path: Path) -> dict[str, str]:
    """<rank>.batch_sample_status.tsv -> {sample: 'included' | 'excluded: <reason>'}.

    `unused` rows (metadata without a count-table column) are not samples of the analysis.
    """
    statuses, seen = {}, set()
    with path.open() as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            sample = (record.get("sample") or "").strip()
            status = (record.get("status") or "").strip()
            reason = (record.get("reason") or "").strip()
            if not sample or status not in ("included", "excluded", "unused"):
                raise ValueError(f"{path.name}: invalid batch ledger row {sample!r} / {status!r}")
            if status == "excluded" and not reason:
                raise ValueError(f"{path.name}: excluded sample {sample} has no reason")
            # Every row counts, unused ones included: an ID may appear once per ledger.
            if sample in seen:
                raise ValueError(f"{path.name}: duplicate sample {sample}")
            seen.add(sample)
            if status != "unused":
                statuses[sample] = "included" if status == "included" else f"excluded: {reason}"
    return statuses


def sample_status(row: dict) -> str:
    def is_zero(key):
        return str(row.get(key)) in ("0", "0.0")
    if row.get("qc_floor") not in (None, NA, "pass"):
        return "not_assessable"
    if is_zero("reads_after_fastp"):
        return "dropped_zero_after_fastp"
    if is_zero("reads_host_depleted"):
        return "no_reads_after_depletion"
    if row.get("bracken_status") == "no_call":
        return "below_bracken_threshold"
    if row.get("bracken_status") == "called":
        return "profiled"
    return NA


def check_low_as(row: dict) -> None:
    """The rescued counts can only shrink stage by stage, and are 0 when the rescue was off."""
    def number(key):
        value = row.get(key)
        return None if value in (None, "", NA) else int(float(value))
    extracted = number("low_as_reads_extracted")
    if row.get("low_as_mode") == "off" and extracted not in (None, 0):
        raise ValueError(f"{row['sample']}: low_as_mode off but low_as_reads_extracted = {extracted}")
    chain = [("low_as_reads_extracted", extracted),
             ("low_as_reads_after_depletion", number("low_as_reads_after_depletion")),
             ("low_as_reads_classified", number("low_as_reads_classified"))]
    known = [(k, v) for k, v in chain if v is not None]
    for (k1, v1), (k2, v2) in zip(known, known[1:]):
        if v2 > v1:
            raise ValueError(f"{row['sample']}: {k2} ({v2}) exceeds {k1} ({v1})")


# --------------------------------------------------------------------------- build


def sample_of(path: Path, suffix: str) -> str:
    return path.name[: -len(suffix)] if path.name.endswith(suffix) else path.stem


def build(inputs: Path, min_reads: int, stages: list[str] | None = None,
          min_after_depletion: int = 0) -> list[dict]:
    rows: dict[str, dict] = {}

    def row_for(sample: str) -> dict:
        return rows.setdefault(sample, {c: NA for c in COLUMNS} | {"sample": sample})

    for path in sorted(inputs.rglob("*.fastp.json")):
        row_for(sample_of(path, ".fastp.json")).update(parse_fastp(path))

    for path in sorted(inputs.rglob("*.extract_counts.tsv")):
        row_for(sample_of(path, ".extract_counts.tsv")).update(parse_extract_counts(path))

    for path in sorted(inputs.rglob("*.upstream_aligner.tsv")):
        row_for(sample_of(path, ".upstream_aligner.tsv")).update(parse_upstream_aligner(path))

    for path in sorted(inputs.rglob("*.input_provenance.tsv")):
        row_for(sample_of(path, ".input_provenance.tsv")).update(parse_input_provenance(path))

    for path in sorted(inputs.rglob("*.low_as_survival.tsv")):
        row_for(sample_of(path, ".low_as_survival.tsv")).update(parse_low_as_survival(path))

    for path in sorted(inputs.rglob("*.depletion_counts.tsv")):
        row_for(sample_of(path, ".depletion_counts.tsv")).update(parse_depletion(path))

    for path in sorted(inputs.rglob("*.krakenuniq.report.txt")):
        sample = sample_of(path, ".krakenuniq.report.txt")
        reads = parse_krakenuniq_microbial_reads(path)
        row = row_for(sample)
        row.update(parse_krakenuniq_clades(path))
        # The gate uses microbial reads (root minus human and synthetic construct), exactly as
        # Bracken.nf applies it (decision 1, 27 Oct 2026).
        micro = row.get("bracken_microbial_reads")
        row["cleared_consensus_gate"] = (
            NA if micro in (None, NA) else ("true" if micro >= min_reads else "false")
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
        row_for(sample).update(parse_bracken_nonmicrobial(path))

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

    # Decontamination and batch correction also CREATE rows, so a run that starts at either
    # stage still gets one row per analysed sample. A metadata-only decontam row does not.
    decontam_output = None   # samples decontamination passed on (final_status included)
    for path in sorted(inputs.rglob("*sample_validation.tsv")):
        if decontam_output is not None:
            raise ValueError("More than one decontamination sample_validation.tsv reached the QC summary")
        decontam_output = set()
        for sample, roles in parse_sample_validation(path).items():
            in_table = roles.pop("_in_table")
            if roles["decontam_final_status"] == "included":
                decontam_output.add(sample)
            if sample in rows or in_table:
                row_for(sample).update(roles)

    ledgers = {}
    for path in sorted(inputs.rglob("*.batch_sample_status.tsv")):
        rank = path.name[: -len(".batch_sample_status.tsv")]
        if rank not in ("species", "genus") or rank in ledgers:
            raise ValueError(f"Unexpected or repeated batch ledger rank in {path.name}")
        ledgers[rank] = parse_batch_sample_status(path)
    # Every rank is corrected from the same samples, and after decontamination from exactly the
    # samples it passed on. Excluded samples are fine when the ledger gives the reason (e.g. outside
    # the tumor subset, all-zero); a sample that appears or vanishes between stages is not.
    sample_sets = {rank: set(ledger) for rank, ledger in ledgers.items()}
    if len({frozenset(s) for s in sample_sets.values()}) > 1:
        raise ValueError(f"Batch ledgers cover different samples per rank: {sorted(set.symmetric_difference(*sample_sets.values()))[:5]}")
    for rank, samples in sample_sets.items():
        if decontam_output is not None and samples != decontam_output:
            raise ValueError(
                f"{rank} batch ledger does not match decontamination output: "
                f"not from decontamination {sorted(samples - decontam_output)[:5]}, "
                f"missing from batch ledger {sorted(decontam_output - samples)[:5]}")
    for rank, ledger in ledgers.items():
        for sample, status in ledger.items():
            row_for(sample)[f"batch_{rank}_status"] = status

    # Samples the sheet asked for that produced nothing at all (audit C10).
    expected: set[str] = set()
    for path in sorted(inputs.rglob("expected_samples.txt")):
        for sample in (s.strip() for s in path.read_text().splitlines()):
            if sample:
                expected.add(sample)
            if sample and sample not in rows:
                row_for(sample)["sample_status"] = "no_outputs"

    # Each sample's input type (sample<TAB>bam|cram|fastq|host_depleted), so extraction is only
    # expected of BAM/CRAM inputs; and the samples whose HUMAnN profile finished.
    input_types: dict[str, str] = {}
    for path in sorted(inputs.rglob("sample_input_types.tsv")):
        for line in path.read_text().splitlines():
            fields = line.split("\t")
            if len(fields) == 2 and fields[0].strip():
                input_types[fields[0].strip()] = fields[1].strip()
    humann_done: set[str] = set()
    for path in sorted(inputs.rglob("humann3_completed.txt")):
        humann_done |= {s.strip() for s in path.read_text().splitlines() if s.strip()}
    if stages is None:
        stages = infer_stages(rows)

    # derived retention, once every count is in
    for row in rows.values():
        # A sample fastp emptied never reaches depletion: its depletion counts are 0, not unknown.
        if str(row.get("reads_after_fastp")) == "0":
            for column in ("reads_after_hg38", "reads_after_t2t_phix", "reads_host_depleted"):
                if row.get(column) in (None, NA):
                    row[column] = 0
        row["qc_floor"] = qc_floor(row, input_types.get(row["sample"]), min_after_depletion)
        recorded = row.pop("_input_qc_floor", None)
        if recorded is not None and recorded != "pass" and row["qc_floor"] != recorded:
            raise ValueError(f"{row['sample']}: extractReads recorded input_qc_floor {recorded}, "
                             f"but the QC summary computes {row['qc_floor']}")
        if recorded == "pass" and row["qc_floor"] in QC_FLOOR_INPUT:
            raise ValueError(f"{row['sample']}: extractReads recorded input_qc_floor pass, "
                             f"but the QC summary computes {row['qc_floor']}")
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
        check_low_as(row)
        if row.get("low_as_mode") == "off":
            # Nothing was rescued, so nothing can survive: 0, not unknown (audit C13).
            for column in ("low_as_reads_after_depletion", "low_as_reads_classified"):
                if row.get(column) in (None, NA) and row.get("low_as_reads_extracted") == "0":
                    row[column] = 0
        row["low_as_share_of_bracken"] = pct(row.get("low_as_reads_classified"),
                                             row.get("kraken_classified_reads"))
        # Gated stages and the terminal state (audit X05). A sample is per-sample when the sheet
        # named it or a per-sample stage reported it; rows created only by decontamination or batch
        # correction (--start_from decontam / batch_correction) are not.
        per_sample = row["sample"] in expected or any(
            _known(row.get(c)) for c in ("reads_after_fastp", "library_primary_records",
                                         "reads_host_depleted", "bracken_status"))
        row["metaphlan_status"] = gated_status("metaphlan", row, stages, humann_done) if per_sample else NA
        row["humann_status"] = gated_status("humann", row, stages, humann_done) if per_sample else NA
        row["terminal_state"] = terminal_state(row, stages if per_sample else [],
                                               input_types.get(row["sample"]), humann_done, per_sample)
        row.pop("_decontam_reason", None)

    return [rows[s] for s in sorted(rows)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path,
                        help="directory of staged per-sample QC artefacts")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-reads", type=int, default=100000,
                        help="--consensus_min_reads, for cleared_consensus_gate")
    parser.add_argument("--min-reads-after-depletion", type=int, default=0,
                        help="--min_reads_after_depletion, for qc_floor (0 = off)")
    parser.add_argument("--stages", default=None,
                        help="comma-separated per-sample stages this run requested, from: "
                             + ", ".join(STAGE_ORDER) + ". main.nf passes them; without it they "
                             "are inferred from the outputs present (a stage every sample failed "
                             "cannot be inferred)")
    args = parser.parse_args()

    stages = parse_stages(args.stages) if args.stages is not None else None
    rows = build(args.inputs, args.min_reads, stages, args.min_reads_after_depletion)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"QC summary: {len(rows)} sample(s) -> {args.output}")
    # Counts only: which terminal states, how many samples each. Never silent about a failure.
    states: dict[str, int] = {}
    for row in rows:
        states[row["terminal_state"]] = states.get(row["terminal_state"], 0) + 1
    for state, n in sorted(states.items()):
        print(f"  {state}: {n}")
    not_assessable = sum(n for s, n in states.items() if s.startswith("not_assessable_"))
    if not_assessable:
        print(f"WARNING: {not_assessable} sample(s) are NOT ASSESSABLE (terminal_state not_assessable_*): "
              "no unmapped reads to profile, unmapped reads filtered upstream, or too few reads after host "
              "depletion. They are not microbial negatives")
    failed = sum(n for s, n in states.items() if s.startswith("failed_at_"))
    if failed:
        print(f"WARNING: {failed} sample(s) failed at a per-sample stage (terminal_state failed_at_*); "
              "their absence from the cohort tables is a failure, not an absence of microbes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
