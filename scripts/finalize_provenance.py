#!/usr/bin/env python3
"""Finish the run record, and write the run report a person actually reads.

The workflow writes the launch half (commit, parameters, reference digests, requested
stages) before any task runs. This adds what only the end of a run knows -- the
outcome, what each process actually cost, the envs this run used and what they pinned,
and SHA-256 digests of the cohort result tables -- and writes both halves out:

  <tracedir>/provenance/<sessionId>_<runName>.json   one per launch, machine-readable
  <tracedir>/RUN_REPORT.txt                          the same run in reading order, appended

`-resume` reuses the session id, so the run name is part of the key: each launch keeps
its own record (audit C06).

Called from the workflow.onComplete handler in nextflow.config. It lives in Python
rather than inline in the config because classes in lib/ are not visible to the config
parser, and because rendering and parsing belong somewhere the test suite can reach.

Ported from pksProfiler (F16). Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

SIZE_UNITS = {"B": 1.0, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}


# --------------------------------------------------------------------------- parsing


def to_bytes(text: str) -> float | None:
    """'1.3 GB' -> 1.3e9. Nextflow writes '-' for a task that never started."""
    if not text or text.strip() in {"-", ""}:
        return None
    match = re.match(r"^\s*([0-9.]+)\s*([KMGT]?B)\s*$", text.strip(), re.IGNORECASE)
    if not match:
        return None
    return float(match.group(1)) * SIZE_UNITS[match.group(2).upper()]


def human_bytes(value: float | None) -> str:
    if value is None:
        return "-"
    for unit in ("TB", "GB", "MB", "KB"):
        if value >= SIZE_UNITS[unit]:
            return f"{value / SIZE_UNITS[unit]:.2f} {unit}"
    return f"{value:.0f} B"


def process_name(task_name: str) -> str:
    """'filterReads (SAMPLE1)' -> 'filterReads'. Drops the per-sample tag."""
    return re.sub(r"\s*\(.*\)\s*$", "", task_name).strip()


def read_trace(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return []          # header-only: the run recorded no completed task
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) == len(header):
            rows.append(dict(zip(header, fields)))
    return rows


def summarise_processes(rows: list[dict]) -> list[dict]:
    """Per-process cost. This is the table that lets resource requests be sized from
    measurement rather than guessed -- see MERGE_CHECKLIST 4.5."""
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[process_name(row.get("name", "?"))].append(row)

    summary = []
    for name, tasks in sorted(grouped.items()):
        peaks = [b for b in (to_bytes(t.get("peak_rss", "")) for t in tasks) if b is not None]
        cpus = {t.get("%cpu", "-") for t in tasks}
        summary.append(
            {
                "process": name,
                "tasks": len(tasks),
                "failed": sum(1 for t in tasks if t.get("status") not in {"COMPLETED", "CACHED"}),
                "peak_rss_max": max(peaks) if peaks else None,
                "peak_rss_median": sorted(peaks)[len(peaks) // 2] if peaks else None,
                "realtime_sample": tasks[0].get("realtime", "-"),
                "pct_cpu_sample": sorted(cpus)[0] if cpus else "-",
            }
        )
    return summary


# Packages named in the human-readable report. The JSON keeps every package; this list
# only decides what RUN_REPORT.txt prints. Widened for D08 (bowtie2, diamond, htslib,
# kraken2, skani, the R stats packages, ConQuR and cqrReg were all missing).
INTERESTING = {
    "samtools", "htslib", "fastp", "fastqc", "minimap2", "bowtie2", "diamond", "multiqc",
    "krakenuniq", "kraken2", "bracken", "metaphlan", "humann", "skani", "r-base",
    "bioconductor-decontam", "bioconductor-phyloseq", "r-vegan", "r-zcompositions",
    "r-quantreg", "r-glmnet", "python", "pandas", "ConQuR", "cqrReg",
}

# R packages installed into the batch-correction env from GitHub/CRAN after the conda
# solve, so conda-meta cannot see them (audit C06/D03/D08).
R_PACKAGES_OF_RECORD = ("ConQuR", "cqrReg")
DESCRIPTION_FIELDS = ("Package", "Version", "RemoteType", "RemoteRepo", "RemoteUsername",
                      "RemoteRef", "RemoteSha", "Packaged")

PACKAGE_FILE = re.compile(r"^(?P<name>.+)-(?P<version>[^-]+)-(?P<build>[^-]+)\.(?:conda|tar\.bz2)$")


def sha256_file(path: Path) -> str:
    """Content SHA-256, streamed. Only the digest leaves this function."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def explicit_digest(lines: list[str]) -> str:
    """SHA-256 of the sorted `url#md5` lines of an explicit spec.

    Computed the same way from a lockfile and from a built env's conda-meta, so an env
    built exactly from its lockfile carries the same value as the lockfile (audit D08).
    """
    body = "\n".join(sorted(set(line.strip() for line in lines if line.strip())))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def lockfile_packages(lock: Path) -> tuple[list[str], list[str]]:
    """(name=version=build list, url#md5 lines) of a `conda list --explicit --md5` lockfile."""
    packages, explicit = [], []
    for line in lock.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("@"):
            continue
        explicit.append(line)
        filename = line.split("#", 1)[0].rsplit("/", 1)[-1]
        match = PACKAGE_FILE.match(filename)
        if match:
            packages.append(f"{match['name']}={match['version']}={match['build']}")
    return sorted(packages), explicit


def read_description(path: Path) -> dict[str, str]:
    """The named fields of an R package DESCRIPTION (Debian control format).

    Package metadata only: version, where it was installed from, and the git commit.
    """
    fields: dict[str, str] = {}
    if not path.is_file():
        return fields
    key = None
    for line in path.read_text(errors="replace").splitlines():
        if line[:1] in (" ", "\t") and key:
            fields[key] += " " + line.strip()
            continue
        if ":" in line:
            key, value = line.split(":", 1)
            key = key.strip()
            fields[key] = value.strip()
    return {k: fields[k] for k in DESCRIPTION_FIELDS if k in fields}


def prefix_env(prefix: Path, conqur_sha: str | None = None) -> dict:
    """A pre-built env addressed by its absolute prefix (params.batch_corr_env).

    Reads conda-meta for what conda installed, and lib/R/library/*/DESCRIPTION for the R
    packages installed after the solve, which conda-meta cannot see. One directory level
    each; nothing else in the prefix is walked.
    """
    entry: dict = {"kind": "prefix", "exists": prefix.is_dir()}
    if not prefix.is_dir():
        return entry
    packages, explicit, owned = [], [], set()
    meta = prefix / "conda-meta"
    for record_file in sorted(meta.glob("*.json")) if meta.is_dir() else []:
        try:
            info = json.loads(record_file.read_text())
        except (OSError, ValueError):
            continue
        name, version, build = info.get("name"), info.get("version"), info.get("build")
        if name and version:
            packages.append(f"{name}={version}={build}")
        if info.get("url"):
            explicit.append(f"{info['url']}#{info['md5']}" if info.get("md5") else info["url"])
        for owned_file in info.get("files") or []:
            parts = owned_file.split("/")
            if len(parts) > 3 and parts[:3] == ["lib", "R", "library"]:
                owned.add(parts[3])
    entry["packages"] = sorted(packages)
    entry["package_count"] = len(packages)
    entry["explicit_sha256"] = explicit_digest(explicit) if explicit else None

    library = prefix / "lib" / "R" / "library"
    if library.is_dir():
        entry["r_packages"] = {}
        for name in R_PACKAGES_OF_RECORD:
            desc = read_description(library / name / "DESCRIPTION")
            entry["r_packages"][name] = desc or "absent"
        # Everything in the R library that no conda package owns: installed by
        # install_conqur.R or at run time, from CRAN/GitHub (audit D03).
        untracked = {}
        for pkg_dir in sorted(library.iterdir()):
            if pkg_dir.is_dir() and pkg_dir.name not in owned and (pkg_dir / "DESCRIPTION").is_file():
                desc = read_description(pkg_dir / "DESCRIPTION")
                untracked[pkg_dir.name] = "@".join(
                    v for v in (desc.get("Version"), desc.get("RemoteSha")) if v) or "?"
        entry["r_packages_not_from_conda"] = untracked
        conqur = entry["r_packages"].get("ConQuR")
        if conqur_sha and isinstance(conqur, dict):
            entry["conqur_sha_matches_param"] = conqur.get("RemoteSha") == conqur_sha
    return entry


def env_record(spec: str, conqur_sha: str | None = None) -> dict:
    """What one conda env spec pinned, and how to tell whether it changed."""
    path = Path(spec)
    if path.is_dir():
        entry = prefix_env(path, conqur_sha)
    elif path.is_file():
        entry = {"kind": "lockfile" if path.name.endswith(".txt") else "yml",
                 "sha256": sha256_file(path)}
        if entry["kind"] == "lockfile":
            packages, explicit = lockfile_packages(path)
            entry["packages"] = packages
            entry["package_count"] = len(packages)
            entry["explicit_sha256"] = explicit_digest(explicit) if explicit else None
            # conda_envs/lock/<env>.linux-64.txt <- conda_envs/<env>.yml, its source.
            stem = path.name.split(".linux-64")[0]
            source = path.parent.parent / f"{stem}.yml"
            if source.is_file():
                entry["source_yml"] = str(source)
                entry["source_yml_sha256"] = sha256_file(source)
    else:
        entry = {"kind": "absent"}
    entry["spec"] = spec
    return entry


def conda_envs(specs: dict[str, str], conqur_sha: str | None = None) -> dict[str, dict]:
    """One entry per env THIS run used, keyed by its parameter name (audit C06/D08).

    `specs` is the resolved params.*_env map the launch half wrote, already narrowed to
    the stages this launch ran. The conda cache directory is not walked: it also holds
    every stale env from earlier yml versions, and envs there are named by an opaque hash.
    """
    return {name: env_record(spec, conqur_sha) for name, spec in sorted((specs or {}).items()) if spec}


def software_summary(envs: dict[str, dict]) -> dict[str, list[str]]:
    """name -> the INTERESTING packages of each env, for RUN_REPORT.txt."""
    summary: dict[str, list[str]] = {}
    for name, entry in envs.items():
        picked = []
        for package in entry.get("packages", []):
            pkg, _, rest = package.partition("=")
            if pkg in INTERESTING:
                picked.append(f"{pkg}={rest.split('=')[0]}")
        for pkg, desc in (entry.get("r_packages") or {}).items():
            if isinstance(desc, dict):
                sha = desc.get("RemoteSha", "")
                picked.append(f"{pkg}={desc.get('Version', '?')}" + (f"@{sha[:8]}" if sha else ""))
        if picked:
            summary[name] = sorted(set(picked))
    return summary


# ------------------------------------------------------------------ output checksums

# Cohort-level result tables, per results-directory parameter: explicit globs, at most
# two levels deep, never a recursive scan of params.outdir (audit C06). Per-sample
# outputs (reads, per-sample Bracken/MetaPhlAn reports, validation per_sample/) are
# deliberately not listed.
COHORT_TABLES = {
    "consensus_taxa_dir": ["*.txt", "*.tsv"],
    "krakenuniq_bracken_dir": ["bracken.genus.mpa.report.txt", "bracken.species.mpa.report.txt"],
    "metaphlan4_dir": ["merged_abundance_table*.txt"],
    "humann3_dir": ["merged/*.tsv"],
    "decontam_dir": ["*/*.csv", "*/*.tsv", "*/*.txt"],
    "batch_corr_dir": ["*/corrected/*", "*/normalized/*.tsv", "*/permanova_summary.tsv",
                       "*/batch_correction_sample_status.tsv"],
    "validation_dir": ["*.tsv"],
    "qc_summary_dir": ["cmp.qc.summary.tsv", "cmp_cohort_report.html"],
}
READ_FILE = re.compile(r"\.(fastq|fq|fasta|fa|fna|sam|bam|cram|bai|crai|bowtie2\.bz2)(\.gz|\.bz2)?$",
                       re.IGNORECASE)
MAX_OUTPUT_FILES = 2000
MAX_OUTPUT_BYTES = 8 * 1024 ** 3   # per file; a cohort table is megabytes


def output_checksums(params: dict) -> dict:
    """SHA-256, size and path relative to params.outdir of each cohort table present.

    Hashes only: no content is read into the record or printed. Paths outside outdir
    are given as <param>/<path within that directory>.
    """
    outdir = Path(params["outdir"]).resolve() if params.get("outdir") else None
    files: list[dict] = []
    seen: set[Path] = set()
    truncated = False
    for param, patterns in COHORT_TABLES.items():
        root_text = params.get(param)
        if not root_text or root_text in {"null", "None"}:
            continue
        root = Path(root_text)
        if not root.is_dir():
            continue
        for pattern in patterns:
            for path in sorted(root.glob(pattern)):
                if not path.is_file() or READ_FILE.search(path.name):
                    continue
                resolved = path.resolve()
                if resolved in seen:
                    continue
                if len(files) >= MAX_OUTPUT_FILES:
                    truncated = True
                    break
                seen.add(resolved)
                try:
                    relative = str(resolved.relative_to(outdir)) if outdir else None
                except ValueError:
                    relative = None
                relative = relative or f"{param}/{path.relative_to(root)}"
                size = path.stat().st_size
                files.append({
                    "path": relative,
                    "bytes": size,
                    "sha256": sha256_file(path) if size <= MAX_OUTPUT_BYTES else "skipped: too large",
                })
    files.sort(key=lambda entry: entry["path"])
    return {"files": files, "file_count": len(files), "truncated": truncated,
            "total_bytes": sum(entry["bytes"] for entry in files)}


# ------------------------------------------------------------------------- rendering


def render(record: dict, processes: list[dict], software: dict[str, list[str]],
           tracedir: Path, session: str) -> str:
    run = record.get("run", {})
    code = record.get("code", {})
    invocation = record.get("invocation", {})
    stats = run.get("task_counts", {})

    out: list[str] = []
    add = out.append
    bar = "=" * 78

    add(bar)
    add(f"CMPipeline run report -- {record.get('status', '?').upper()}")
    add(bar)
    add(f"  session       {session}")
    add(f"  started       {run.get('started', '-')}")
    add(f"  completed     {run.get('completed', '-')}")
    add(f"  duration      {run.get('duration', '-')}")
    add(f"  exit status   {run.get('exit_status', '-')}")
    if record.get("status") == "failed" and run.get("error_message"):
        add(f"  error         {run['error_message'].splitlines()[0][:70]}")
    add("")

    add("REPRODUCIBILITY")
    add(f"  commit        {code.get('commit', '-')} ({code.get('commit_source', '-')})")
    if code.get("dirty"):
        add("  *** UNCOMMITTED CHANGES -- this run is NOT attributable to that commit ***")
        add(f"      {len(code.get('modified_files', []))} modified path(s) in the working tree")
    else:
        add("  working tree  clean")
    add(f"  scripts       {code.get('scripts_digest', '-')}")
    add(f"  modules       {code.get('modules_digest', '-')}")
    add(f"  conda_envs/   {code.get('envs_digest', '-')}")
    add(f"  nextflow      {record.get('environment', {}).get('nextflow_version', '-')}")
    add("")

    add("WHAT WAS REQUESTED")
    for stage in invocation.get("stages_requested", []):
        add(f"  - {stage}")
    add("")

    add("WHAT RAN")
    if stats:
        add(f"  succeeded {stats.get('succeeded', '-')}   cached {stats.get('cached', '-')}"
            f"   failed {stats.get('failed', '-')}   ignored {stats.get('ignored', '-')}")
    if processes:
        add("")
        add(f"  {'process':<22}{'tasks':>6}{'failed':>7}{'peak RSS max':>15}{'median':>12}")
        for entry in processes:
            add(f"  {entry['process']:<22}{entry['tasks']:>6}{entry['failed']:>7}"
                f"{human_bytes(entry['peak_rss_max']):>15}"
                f"{human_bytes(entry['peak_rss_median']):>12}")
        add("")
        add("  Peak RSS is what resource requests should be sized from. If a request in")
        add("  conf/base.config is orders of magnitude above the peak here, it is costing")
        add("  throughput on a shared cluster.")
    else:
        add("  per-process detail unavailable: pipeline_info/trace.txt had no task rows")
    add("")

    failures = [entry for entry in processes if entry["failed"]]
    if failures:
        add("WHAT FAILED")
        for entry in failures:
            add(f"  {entry['process']}: {entry['failed']} of {entry['tasks']} task(s)")
        add("")

    add("REFERENCES AND DATABASES")
    for name, digest in (record.get("dependencies") or {}).items():
        add(f"  {name:<22}{digest}")
    add("")

    envs = (record.get("environment") or {}).get("conda_envs_used") or {}
    if envs:
        add("CONDA ENVS THIS RUN USED")
        for env_name, entry in envs.items():
            digest = entry.get("explicit_sha256") or entry.get("sha256") or "-"
            add(f"  {env_name:<24}{entry.get('kind', '?'):<10}{str(digest)[:12]}  {entry.get('spec', '-')}")
            if software.get(env_name):
                add(f"      {', '.join(software[env_name])}")
            extra = entry.get("r_packages_not_from_conda")
            if extra:
                add(f"      {len(extra)} R package(s) not installed by conda")
            if entry.get("conqur_sha_matches_param") is False:
                add("      *** ConQuR RemoteSha differs from --conqur_sha ***")
        add("")

    tables = (record.get("outputs") or {}).get("cohort_tables") or {}
    add("COHORT TABLES (sha256)")
    if tables.get("files"):
        for entry in tables["files"]:
            add(f"  {str(entry['sha256'])[:16]}  {human_bytes(entry['bytes']):>10}  {entry['path']}")
        if tables.get("truncated"):
            add(f"  ... truncated at {tables['file_count']} files")
    else:
        add("  none found under the results directories")
    add("")

    add("WHERE")
    add(f"  {'full record':<22}{tracedir}/provenance/{session}.json")
    add(f"  {'trace':<22}{tracedir}/trace.txt")
    add(f"  {'report / timeline':<22}{tracedir}/report.html, {tracedir}/timeline.html")
    add(f"  {'work dir':<22}{invocation.get('work_dir', '-')}")
    add(bar)
    add("")
    return "\n".join(out)


# ------------------------------------------------------------------------------ main


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tracedir", required=True, type=Path)
    parser.add_argument("--session", required=True)
    parser.add_argument("--success", default="false")
    parser.add_argument("--exit-status", default="")
    parser.add_argument("--duration", default="")
    parser.add_argument("--error-message", default="")
    parser.add_argument("--succeeded", default="0")
    parser.add_argument("--cached", default="0")
    parser.add_argument("--failed", default="0")
    parser.add_argument("--ignored", default="0")
    # Accepted and ignored: the cache dir is no longer walked (audit C06/D08).
    parser.add_argument("--conda-cache-dir", default="", help=argparse.SUPPRESS)
    args = parser.parse_args()

    target = args.tracedir / "provenance" / f"{args.session}.json"
    if not target.is_file():
        # The launch half never got written; say so rather than inventing a record.
        print(f"WARN: no launch record at {target}; nothing to finalise")
        return 0
    record = json.loads(target.read_text())

    success = args.success.strip().lower() == "true"
    record["status"] = "completed" if success else "failed"
    record.setdefault("run", {}).update(
        {
            "completed": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "duration": args.duration,
            "exit_status": args.exit_status,
            "error_message": args.error_message,
            "task_counts": {
                "succeeded": args.succeeded,
                "cached": args.cached,
                "failed": args.failed,
                "ignored": args.ignored,
            },
        }
    )

    rows = read_trace(args.tracedir / "trace.txt")
    processes = summarise_processes(rows)
    launch_params = (record.get("invocation") or {}).get("params") or {}
    record["outputs"] = {
        "processes": processes,
        # Digests of the cohort tables as they stand at the end of this launch (C06).
        "cohort_tables": output_checksums(launch_params),
    }

    environment = record.setdefault("environment", {})
    envs = conda_envs(environment.get("conda_envs") or {}, launch_params.get("conqur_sha"))
    environment["conda_envs_used"] = envs
    software = software_summary(envs)
    if software:
        environment["software"] = software

    target.write_text(json.dumps(record, indent=2) + "\n")

    report = render(record, processes, software, args.tracedir, args.session)
    report_path = args.tracedir / "RUN_REPORT.txt"
    with report_path.open("a") as handle:
        handle.write(report)
    print(f"Run report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
