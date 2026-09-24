#!/usr/bin/env python3
"""Assertions for the executed integration test (audit C22). Test use only.

Reads ONLY the synthetic run's outputs under the fixture directory and the expected counts
make_fixtures.py wrote. Prints one PASS/FAIL line per assertion, with counts only, and exits
with the number of failures (0 = all passed).

    python3 check_outputs.py --fixtures <dir> --run-a <dir> --run-b <dir> --batch-correction on|off
"""

import argparse
import csv
import json
import re
from pathlib import Path

REQUIRED_QC_COLUMNS = [
    "sample", "library_primary_records", "extracted_unmapped_records", "pct_library_unmapped",
    "reads_raw", "reads_after_fastp", "reads_after_hg38", "reads_after_t2t_phix",
    "reads_host_depleted", "pct_retained_overall", "pct_retained_of_library",
    "bracken_microbial_reads", "cleared_consensus_gate", "bracken_genera", "bracken_species",
    "decontam_sample_type", "decontam_role", "decontam_final_status",
]

results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def exit_status(run):
    try:
        return int((run / "exit_status").read_text().strip())
    except (OSError, ValueError):
        return None


def task_dirs(run, process):
    """Work dirs of every task of `process` in this run, from pipeline_info/trace.txt."""
    trace = run / "pipeline_info/trace.txt"
    if not trace.exists():
        return []
    dirs = []
    for row in read_tsv(trace):
        if re.match(rf"^{re.escape(process)}( \(|$)", row.get("name", "")):
            dirs += sorted((run / "work").glob(row["hash"] + "*"))
    return dirs


def log_has(run, text):
    for p in (run / ".nextflow.log", run / "nextflow.stdout"):
        if p.exists() and text in p.read_text(errors="replace"):
            return True
    return False


def run_a(fx, run, expected):
    res = run / "results"
    samples = expected["samples"]
    check("A: pipeline exit status 0", exit_status(run) == 0, f"exit={exit_status(run)}")
    if not (run / "pipeline_info/trace.txt").exists():
        # Nothing ran (e.g. a compile error): one clear failure, not 70 derivative ones.
        check("A: pipeline started any task (see runA/nextflow.stdout)", False, "no trace.txt")
        return

    # ---- extraction: library denominator and every unmapped category (audit C07, F01)
    for s, e in samples.items():
        if e["kind"] not in ("bam", "cram"):
            continue
        f = res / f"UNMAPPED_BAM/{s}.extract_counts.tsv"
        rows = read_tsv(f) if f.exists() else []
        check(f"A: {e['kind']} extract_counts.tsv exists", len(rows) == 1)
        if rows:
            got_p, got_u = rows[0].get("library_primary_records"), rows[0].get("extracted_unmapped_records")
            check(f"A: {e['kind']} library_primary_records == {e['library_primary_records']}",
                  got_p == str(e["library_primary_records"]), f"got {got_p}")
            check(f"A: {e['kind']} extracted_unmapped_records == {e['extracted_unmapped_records']}",
                  got_u == str(e["extracted_unmapped_records"]), f"got {got_u}")
    leaked = [p.name for p in (res / "UNMAPPED_BAM").glob("*.f*q.gz")] if (res / "UNMAPPED_BAM").exists() else []
    check("A: no intermediate read files published without --save_intermediates (C11)",
          not leaked, f"{len(leaked)} read files")

    # ---- host depletion counts, stage by stage
    for s, e in samples.items():
        f = res / f"MAPPED_READS/{s}.depletion_counts.tsv"
        stages = {r["stage"]: r["reads"] for r in read_tsv(f)} if f.exists() else {}
        want = {"input": e["reads_raw"], "after_hg38": e["reads_after_hg38"],
                "after_t2t_phix": e["reads_after_t2t_phix"], "host_depleted": e["reads_host_depleted"]}
        got = {k: stages.get(k) for k in want}
        check(f"A: {s} depletion counts {'/'.join(str(v) for v in want.values())}",
              got == {k: str(v) for k, v in want.items()},
              "got " + "/".join(str(v) for v in got.values()))
        check(f"A: {s} host_depleted.fastq.gz published", (res / f"MAPPED_READS/{s}.host_depleted.fastq.gz").exists())

    # ---- Bracken: real Bracken on the stub report; threshold is a param, not task.cpus (C01)
    for s in samples:
        check(f"A: {s} KrakenUniq report + Bracken G/S reports",
              all((res / f"BRACKEN/{s}.{x}").exists() for x in
                  ("krakenuniq.report.txt", "bracken.G.report.txt", "bracken.S.report.txt",
                   "bracken.G.mpa.krakenreport.txt", "bracken.S.mpa.krakenreport.txt")))
    thin = [s for s, e in samples.items() if e["microbial_reads"] < 2]
    for s in thin:
        f = res / f"BRACKEN/{s}.bracken.G.report.txt"
        n = len(f.read_text().splitlines()) if f.exists() else -1
        check("A: thin sample Bracken below-threshold no-call is header-only (C02)", n == 1, f"lines={n}")
    scripts = [d / ".command.sh" for d in task_dirs(run, "Bracken")]
    rendered = [p.read_text() for p in scripts if p.exists()]
    thr = [set(re.findall(r"bracken .*? -t (\S+)", t)) for t in rendered]
    check("A: Bracken rendered with -t 2 (--bracken_threshold), not the CPU count (C01)",
          rendered and all(t == {"2"} for t in thr), f"{len(rendered)} tasks, values {sorted(set().union(*thr)) if thr else []}")
    check("A: merged Bracken genus/species tables",
          (res / "BRACKEN/bracken.genus.mpa.report.txt").exists() and (res / "BRACKEN/bracken.species.mpa.report.txt").exists())

    # ---- MetaPhlAn gate and the consensus pass-through
    check("A: no MetaPhlAn task ran (every sample below --consensus_min_reads)", not task_dirs(run, "metaphlan4"))
    st = res / "CONSENSUS_TAXA/consensus_status.tsv"
    status = {r["key"]: r["value"] for r in read_tsv(st)} if st.exists() else {}
    check("A: consensus_status mode == bracken_passthrough", status.get("mode") == "bracken_passthrough",
          f"mode={status.get('mode')}")
    check(f"A: consensus saw {len(samples)} Bracken samples",
          status.get("samples_in_bracken") == str(len(samples)), f"got {status.get('samples_in_bracken')}")
    lib = res / "CONSENSUS_TAXA/bracken_library_sizes.tsv"
    libs = {r["sample"]: r for r in read_tsv(lib)} if lib.exists() else {}
    nonthin = [s for s in samples if s not in thin]
    check("A: bracken_library_sizes.tsv > 0 for every non-thin sample",
          all(float(libs.get(s, {}).get("bracken_genus_total") or 0) > 0 for s in nonthin),
          f"{len(libs)} rows")
    check("A: consensus genus table written", (res / "CONSENSUS_TAXA/bracken.metaphlan.common.genus.mpa.report.txt").exists())

    # ---- decontamination on the consensus table, wired with its library sizes
    d = res / "04_DECONTAMINATION/consensus_run"
    for suffix in ("contaminants_final.tsv", "sample_validation.tsv", "decontam_pkg_decontaminated.csv",
                   "batch_eligibility.tsv"):
        check(f"A: decontam consensus_run.{suffix}", (d / f"consensus_run.{suffix}").exists())
    elig = read_tsv(d / "consensus_run.batch_eligibility.tsv") if (d / "consensus_run.batch_eligibility.tsv").exists() else []
    check("A: decontam saw 2 eligible batches", sum(r.get("eligible") == "TRUE" for r in elig) == 2, f"{len(elig)} batches")

    # ---- QC summary: column contract and values
    qc = res / "QC_SUMMARY/cmp.qc.summary.tsv"
    rows = read_tsv(qc) if qc.exists() else []
    cols = list(rows[0].keys()) if rows else []
    missing = [c for c in REQUIRED_QC_COLUMNS if c not in cols]
    check("A: QC summary has every required column", qc.exists() and not missing, f"missing={missing}")
    check(f"A: QC summary has {len(samples)} rows", len(rows) == len(samples), f"rows={len(rows)}")
    by = {r["sample"]: r for r in rows}
    for s, e in samples.items():
        r = by.get(s, {})
        for col in ("library_primary_records", "extracted_unmapped_records", "reads_raw", "reads_host_depleted"):
            check(f"A: QC {s} {col} == {e[col]}", r.get(col) == str(e[col]), f"got {r.get(col)}")
        check(f"A: QC {s} bracken_microbial_reads == {e['microbial_reads']}",
              r.get("bracken_microbial_reads") == str(e["microbial_reads"]), f"got {r.get('bracken_microbial_reads')}")
        check(f"A: QC {s} cleared_consensus_gate == false", r.get("cleared_consensus_gate") == "false",
              f"got {r.get('cleared_consensus_gate')}")
        check(f"A: QC {s} decontam_role filled", r.get("decontam_final_status") not in (None, "", "NA"),
              f"status={r.get('decontam_final_status')}")

    check("A: MultiQC report", (res / "MULTIQC/multiqc_report.html").exists())
    check("A: cohort report", (res / "QC_SUMMARY/cmp_cohort_report.html").exists())
    prov = list((run / "pipeline_info/provenance").glob("*.json")) if (run / "pipeline_info/provenance").exists() else []
    check("A: provenance record written", len(prov) >= 1, f"{len(prov)} files")


def run_b(fx, run, expected, batch_correction):
    res = run / "results"
    e = expected["decontam_run_b"]
    code = exit_status(run)
    if not (run / "pipeline_info/trace.txt").exists():
        check("B: exit status 0" if batch_correction == "off" else "B: pipeline ran", False, f"exit={code}, no trace.txt")
        return
    if batch_correction == "on":
        # 12 samples < the 20 run_auto_phase needs: ConQuR must refuse loudly, never pass raw
        # counts off as corrected (audit A01).
        check("B: exit status non-zero (batch correction must fail on 12 samples)", code not in (0, None), f"exit={code}")
        check("B: failure names 'Insufficient samples'", log_has(run, "Insufficient samples"))
        check("B: no ConQuR_tuned.tsv published", not list((res / "05_BATCH_CORRECTION").rglob("ConQuR_tuned.tsv"))
              if (res / "05_BATCH_CORRECTION").exists() else True)
    else:
        check("B: exit status 0", code == 0, f"exit={code}")
    d = res / "04_DECONTAMINATION/decontam_run"
    for suffix in ("contaminants_final.tsv", "sample_validation.tsv", "decontam_pkg_decontaminated.csv",
                   "batch_eligibility.tsv", "validated_metadata.tsv"):
        check(f"B: decontam decontam_run.{suffix}", (d / f"decontam_run.{suffix}").exists())
    final = read_tsv(d / "decontam_run.contaminants_final.tsv") if (d / "decontam_run.contaminants_final.tsv").exists() else []
    check(f"B: contaminants_final evaluates {e['taxa']} taxa", len(final) == e["taxa"], f"rows={len(final)}")
    sv = read_tsv(d / "decontam_run.sample_validation.tsv") if (d / "decontam_run.sample_validation.tsv").exists() else []
    check(f"B: all {e['samples']} samples included", sum(r.get("final_status") == "included" for r in sv) == e["samples"],
          f"included={sum(r.get('final_status') == 'included' for r in sv)}")
    flagged = sum(r.get("final_contaminant") == "TRUE" for r in final)
    print(f"INFO  B: taxa flagged as contaminant: {flagged} of {len(final)}")


def run_c(fx, run, expected, batch_correction):
    """Audit C22: batch correction on enough samples must produce a real correction."""
    if batch_correction != "on":
        print("SKIP  C: batch-correction env not available")
        return
    check("C: exit status 0", exit_status(run) == 0, f"exit={exit_status(run)}")
    status = sorted(run.glob("results/05_BATCH_CORRECTION/**/BATCH_CORRECTION_STATUS.tsv"))
    rows = read_tsv(status[0]) if status else []
    check("C: BATCH_CORRECTION_STATUS says corrected", bool(rows) and rows[0].get("status") == "corrected",
          f"status={rows[0].get('status') if rows else None}")
    check("C: ConQuR_tuned.tsv written", bool(sorted(run.glob("results/05_BATCH_CORRECTION/**/ConQuR_tuned.tsv"))))
    check("C: no UNCORRECTED_passthrough.tsv", not sorted(run.glob("results/**/UNCORRECTED_passthrough.tsv")))


def run_d(fx, run):
    """Audit C04: a missing covariate stops at launch and is named; no task runs."""
    check("D: launch-time stop (non-zero exit)", exit_status(run) not in (None, 0), f"exit={exit_status(run)}")
    check("D: error names the missing covariate", log_has(run, "lacks: bmi"))
    check("D: no task ran", not (run / "pipeline_info/trace.txt").exists() or len(read_tsv(run / "pipeline_info/trace.txt")) == 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", required=True)
    ap.add_argument("--run-a", required=True)
    ap.add_argument("--run-b", required=True)
    ap.add_argument("--batch-correction", choices=["on", "off"], default="off")
    ap.add_argument("--run-c", default=None)
    ap.add_argument("--run-d", default=None)
    args = ap.parse_args()
    expected = json.loads((Path(args.fixtures) / "expected.json").read_text())
    run_a(Path(args.fixtures), Path(args.run_a), expected)
    run_b(Path(args.fixtures), Path(args.run_b), expected, args.batch_correction)
    if args.run_c:
        run_c(Path(args.fixtures), Path(args.run_c), expected, args.batch_correction)
    if args.run_d:
        run_d(Path(args.fixtures), Path(args.run_d))
    failed = results.count(False)
    print(f"SUMMARY  {len(results) - failed} passed, {failed} failed, {len(results)} assertions")
    raise SystemExit(min(failed, 255))


if __name__ == "__main__":
    main()
