#!/usr/bin/env python3
"""Finish the run record, and write the run report a person actually reads.

The workflow writes the launch half (commit, parameters, reference digests, requested
stages) before any task runs. This adds what only the end of a run knows -- the
outcome, what each process actually cost, the envs this run used and what they pinned,
what each task actually ran (the runtime tool-version records, audit K06/K07), and SHA-256
digests of the cohort result tables -- and writes both halves out:

  <tracedir>/provenance/<sessionId>_<runName>.json   one per launch, machine-readable
  <tracedir>/RUN_REPORT.txt                          the same run in reading order, appended

`-resume` reuses the session id, so the run name is part of the key: each launch keeps
its own record (audit C06).

Called from the workflow.onComplete handler in nextflow.config. It lives in Python
rather than inline in the config because classes in lib/ are not visible to the config
parser, and because rendering and parsing belong somewhere the test suite can reach.

Exit status: 0, or 1 when --tool_version_check is strict and a task's record shows a tool that
is not the pinned build (or a recording task left no record). The onComplete handler in main.nf
reports a non-zero exit as an error; it cannot change the run's own status.

Ported from pksProfiler (F16). Standard library only.
"""

from __future__ import annotations

import argparse
import csv
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
    "r-conqur", "r-cqrreg",   # the same two, as the local-channel conda packages (lockfile env)
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


# ------------------------------------------------------------- tool versions (K06/K07)
#
# conda_envs() above is what each env's spec PINNED. What RAN is recorded by the task itself:
# scripts/record_tool_versions.sh, rendered into every process that runs a tool by
# lib/ToolVersions.groovy, asks each binary on the task's PATH -- after conda activation -- for
# its version and location, and writes the answer to the work dir's .tool_versions.tsv. The two
# differ exactly when something earlier on PATH shadows the environment: CRC50 attempt 1 (Job
# 12380606) ran minimap2 2.24 in every mapReads task while the lockfile said 2.28.
# Ported from pksProfiler (F16, 2bc7ca9 / 2c0f8c8), adapted to explicit lockfiles.

TOOL_RECORD = ".tool_versions.tsv"
# Statuses the task itself assigns that mean "not the pinned build from its own env".
BLOCKING_STATUSES = ("OUTSIDE_ENV", "VERSION_MISMATCH", "NOT_FOUND", "NO_ENV")
# A task that finished without a record, although its process writes one.
MISSING = "MISSING"
# The binary that ran reports a version the environment's conda-meta does not hold.
INSTALLED_MISMATCH = "INSTALLED_MISMATCH"
RECORDED_STATUSES = ("COMPLETED", "CACHED")
PROCESS_RE = re.compile(r"^\s*process\s+(\w+)\s*\{", re.M)
ALIAS_RE = re.compile(r"\b(\w+)\s+as\s+(\w+)\b")
ALIASES_FILE = Path(__file__).resolve().parent / "tool_version_aliases.tsv"
META_FILE = re.compile(r"^(?P<name>.+)-(?P<version>[^-]+)-(?P<build>[^-]+)\.json$")


def version_matches(observed: str, expected: str) -> bool:
    """True when `observed` is `expected` or a more specific release of it (x.y.z of x.y)."""
    return observed == expected or observed.startswith(expected + ".")


def version_aliases(path: Path = ALIASES_FILE) -> dict:
    """{(package, conda_version): reported_version} from tool_version_aliases.tsv."""
    aliases = {}
    if not Path(path).is_file():
        return aliases
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) >= 3 and fields[0] and fields[1] and fields[2]:
            aliases[(fields[0], fields[1])] = fields[2]
    return aliases


def version_satisfies(observed: str, expected: str, package: str, aliases: dict) -> bool:
    """version_matches, or the known self-report of `package` at conda version `expected`."""
    if version_matches(observed, expected):
        return True
    alias = aliases.get((package, expected))
    return bool(alias) and version_matches(observed, alias)


def env_pins(spec: str) -> dict[str, str]:
    """{package: exact version} the task's env spec pins, read the way record_tool_versions.sh
    reads it: an explicit lockfile's package URLs, a yml's `- package=version` lines, or an
    existing prefix's conda-meta."""
    pins: dict[str, str] = {}
    path = Path(spec) if spec else None
    if path is None:
        return pins
    if path.is_dir():
        meta = path / "conda-meta"
        for record_file in sorted(meta.glob("*.json")) if meta.is_dir() else []:
            match = META_FILE.match(record_file.name)
            if match:
                pins.setdefault(match["name"], match["version"])
        return pins
    if not path.is_file():
        return pins
    if path.suffix in (".yml", ".yaml"):
        for line in path.read_text().splitlines():
            entry = line.split("#", 1)[0].strip()
            if not entry.startswith("-"):
                continue
            entry = entry[1:].replace(" ", "").split("::")[-1]
            match = re.match(r"^([A-Za-z0-9_.+-]+?)(==?)([^=<>!~,|]+)", entry)
            if match:
                pins.setdefault(match.group(1), match.group(3).rstrip("*").rstrip("."))
        return pins
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("@"):
            continue
        match = PACKAGE_FILE.match(line.split("#", 1)[0].rsplit("/", 1)[-1])
        if match:
            pins.setdefault(match["name"], match["version"])
    return pins


def recording_processes(repo_dir: Path) -> set[str]:
    """Names, as the trace shows them, of the processes that write a tool-version record: every
    process whose body calls ToolVersions.record, and every `include { declared as alias }` of one
    (the trace records the alias, e.g. FASTQC1)."""
    names: set[str] = set()
    sources = sorted(Path(repo_dir).glob("Modules/*.nf")) + [Path(repo_dir) / "main.nf"]
    texts = [p.read_text() for p in sources if p.is_file()]
    for text in texts:
        starts = list(PROCESS_RE.finditer(text))
        for index, match in enumerate(starts):
            end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
            if "ToolVersions.record(" in text[match.end():end]:
                names.add(match.group(1))
    for text in texts:
        for include in re.findall(r"include\s*\{([^}]*)\}", text):
            for declared, alias in ALIAS_RE.findall(include):
                if declared in names:
                    names.add(alias)
    return names


def task_work_dir(row: dict, work_dir: str | None) -> Path | None:
    """The task's work directory: the trace's workdir column, else found from its hash."""
    workdir = (row.get("workdir") or "").strip()
    if workdir and workdir != "-":
        return Path(workdir)
    task_hash = (row.get("hash") or "").strip()
    if work_dir and "/" in task_hash:
        hits = sorted(Path(work_dir).glob(f"{task_hash}*"))
        if len(hits) == 1:
            return hits[0]
    return None


def installed_version(conda_prefix: str, package: str, inventory: dict) -> str | None:
    """What conda-meta says the task's own environment holds for `package`."""
    if not conda_prefix:
        return None
    if conda_prefix not in inventory:
        inventory[conda_prefix] = env_pins(conda_prefix) if Path(conda_prefix).is_dir() else {}
    return inventory[conda_prefix].get(package)


def collect_tool_versions(rows: list[dict], work_dir: str | None, repo_dir: Path,
                          mode: str = "strict") -> dict:
    """Read every task's record, compare it with its pin and with what conda installed, and say
    what does not agree.

    Returns {"check", "tools", "problems", "tasks_recorded", "tasks_unreadable"}. `tools` has one
    row per distinct (tool, version, path, environment); `problems` one per distinct problem, with
    the processes and task count it covers. Tool names, versions and paths only: no sample ID.
    """
    summary = {"check": mode, "tools": [], "problems": [],
               "tasks_recorded": 0, "tasks_unreadable": 0}
    recording = recording_processes(repo_dir)
    aliases = version_aliases()
    pins_by_env: dict[str, dict] = {}
    inventory: dict[str, dict] = {}
    tools: dict[tuple, dict] = {}
    problems: dict[tuple, dict] = {}

    def problem(status, key, detail, process, waived=False):
        # waived: the task ran with a per-task `ext.tool_version_check = 'warn'` (its record's
        # `check` column). Listed in RUN_REPORT, but it does not fail a strict run on its own.
        entry = problems.setdefault((status,) + key, {"status": status, "detail": detail,
                                                      "processes": set(), "tasks": 0, "waived": True})
        entry["processes"].add(process)
        entry["tasks"] += 1
        entry["waived"] = entry["waived"] and waived

    for row in rows:
        process = process_name(row.get("name", "")).rsplit(":", 1)[-1]
        status = (row.get("status") or "").strip()
        directory = task_work_dir(row, work_dir)
        record = directory / TOOL_RECORD if directory else None
        if record is None or not directory.is_dir():
            if process in recording and status in RECORDED_STATUSES:
                summary["tasks_unreadable"] += 1
            continue
        if not record.is_file():
            if process in recording and status in RECORDED_STATUSES:
                problem(MISSING, (process,),
                        f"{process}: finished without a tool-version record", process)
            continue
        summary["tasks_recorded"] += 1
        with record.open(newline="") as handle:
            for entry in csv.DictReader(handle, delimiter="\t"):
                tool = entry.get("tool", "")
                package = entry.get("package") or tool
                version = entry.get("version", "")
                path = entry.get("path", "")
                prefix = entry.get("conda_prefix", "")
                env_file = entry.get("env", "")
                if env_file not in pins_by_env:
                    pins_by_env[env_file] = env_pins(env_file)
                pinned = pins_by_env[env_file].get(package, entry.get("pinned", ""))
                installed = installed_version(prefix, package, inventory) or ""
                task_status = entry.get("status", "")

                # The task judged against the spec it was given; judge again here, against the
                # same spec, and against what conda installed in the prefix that ran.
                verdict = task_status
                if task_status not in BLOCKING_STATUSES:
                    if pinned and version and not version_satisfies(version, pinned, package, aliases):
                        verdict = "VERSION_MISMATCH"
                    elif installed and version and not (
                            version_matches(installed, version)
                            or version_satisfies(version, installed, package, aliases)):
                        verdict = INSTALLED_MISMATCH

                key = (tool, version, path, prefix)
                out = tools.setdefault(key, {
                    "tool": tool, "package": package, "version": version, "pinned": pinned,
                    "installed": installed, "path": path, "conda_prefix": prefix, "env": env_file,
                    "status": verdict, "raw": entry.get("raw", ""), "processes": set(), "tasks": 0})
                out["processes"].add(process)
                out["tasks"] += 1

                if verdict in BLOCKING_STATUSES:
                    where = f"from {path}" if path else "not on PATH"
                    detail = {
                        "OUTSIDE_ENV": f"{tool} {version or '(version unread)'} ran {where}, "
                                       f"outside the task's conda env {prefix}"
                                       + (f"; pinned {package}={pinned}" if pinned else ""),
                        "VERSION_MISMATCH": f"{tool} reports {version} {where}; "
                                            f"{Path(env_file).name or 'its env'} pins {package}={pinned}",
                        "NOT_FOUND": f"{tool} was not on the task's PATH",
                        "NO_ENV": f"{tool} {version} ran {where} with no conda env active",
                    }[verdict]
                    problem(verdict, key, detail, process, waived=entry.get("check") == "warn")
                elif verdict == INSTALLED_MISMATCH:
                    problem(verdict, key, f"{tool} reports {version} from {path}, but conda "
                                          f"installed {package} {installed} in {prefix}", process,
                            waived=entry.get("check") == "warn")

    summary["tools"] = [dict(row, processes=sorted(row["processes"]))
                        for _key, row in sorted(tools.items())]
    summary["problems"] = [dict(entry, processes=sorted(entry["processes"]))
                           for _key, entry in sorted(problems.items())]
    return summary


def render_tool_versions(summary: dict) -> list[str]:
    """RUN_REPORT's TOOL VERSIONS (as run) section, and TOOL VERSION MISMATCH when there is one."""
    out = ["TOOL VERSIONS (as run)"]
    tools = summary.get("tools") or []
    if not tools:
        out.append("  no task left a tool-version record (a run from before the record existed, "
                   "or no trace)")
    else:
        out.append("  reported by each binary from inside its task, after conda activation"
                   " (per task: <work dir>/.tool_versions.tsv)")
        out.append(f"    {'tool':<14}{'version':<11}{'pinned':<11}{'installed':<11}{'tasks':>6}   ran from")
        for row in tools:
            prefix = (row.get("conda_prefix") or "").rstrip("/")
            path = row.get("path") or "-"
            where = f"env {Path(prefix).name}" if prefix and path.startswith(prefix + "/") else path
            flag = "" if row["status"] in ("ok", "unpinned", "version_unread") else f"   {row['status']}"
            out.append(f"    {row['tool']:<14}{row['version'] or '?':<11}{row['pinned'] or '-':<11}"
                       f"{row['installed'] or '-':<11}{row['tasks']:>6}   {where}{flag}")
    if summary.get("tasks_unreadable"):
        out.append(f"  {summary['tasks_unreadable']} task(s) could not be checked: their work"
                   " directory was not found (cleaned up, or another machine)")
    out.append("")

    problems = summary.get("problems") or []
    if problems:
        out.append("TOOL VERSION MISMATCH")
        if summary.get("check") == "strict":
            out.append("  --tool_version_check strict: a task that finds one of these stops the run"
                       " (exit 78); MISSING / INSTALLED_MISMATCH are found here, after the run")
        else:
            out.append("  --tool_version_check warn: the run carried on; these results were NOT"
                       " produced by the pinned tools")
        for entry in problems:
            out.append(f"  {entry['status']:<20}{entry['detail']}")
            out.append(f"  {'':<20}{', '.join(entry['processes'])}; {entry['tasks']} task(s)"
                       + ("; task ran with ext.tool_version_check = 'warn' (listed, not blocking)"
                          if entry.get("waived") else ""))
        out.append("")
    return out


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
                       "*/batch_correction_sample_status.tsv", "*/batch_correction_taxon_status.tsv",
                       "*/batch_correction_quality.tsv", "*/batch_correction_methods.tsv",
                       "*/batch_correction_run_settings.tsv"],
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


# ------------------------------------------------------------- low-score rescue (L8)

LOW_AS_COUNTS = ("low_as_reads_extracted", "low_as_reads_after_depletion", "low_as_reads_classified")


def low_score_rescue_summary(params: dict) -> dict | None:
    """Cohort view of --extract_low_score_mapped from the QC summary: how many BAM/CRAM inputs had
    the rescue on and why, which upstream aligners, and the rescued reads summed per stage.

    Counts only -- no sample ID, path or read enters the record (per-sample values stay in
    cmp.qc.summary.tsv on the cluster). None when there is no QC summary or no BAM/CRAM row.
    """
    root = params.get("qc_summary_dir")
    if not root or root in {"null", "None"}:
        return None
    path = Path(root) / "cmp.qc.summary.tsv"
    if not path.is_file():
        return None
    with path.open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r.get("low_as_mode") in ("on", "off")]
    if not rows:
        return None

    def tally(key):
        out: dict[str, int] = defaultdict(int)
        for r in rows:
            out[key(r)] += 1
        return dict(sorted(out.items()))

    totals = {}
    for column in LOW_AS_COUNTS:
        values = [r.get(column) for r in rows if r.get("low_as_mode") == "on"]
        numbers = [int(float(v)) for v in values if v not in (None, "", "NA")]
        totals[column] = sum(numbers) if numbers else None
        totals[f"{column}_samples"] = len(numbers)
    return {
        "inputs": len(rows),
        "by_mode": tally(lambda r: r["low_as_mode"]),
        "by_reason": tally(lambda r: r.get("low_as_mode_reason") or "NA"),
        "by_aligner": tally(lambda r: f"{r.get('upstream_aligner') or 'NA'}:{r.get('upstream_aligner_version') or 'NA'}"
                                      f" -T {r.get('upstream_bwa_min_score') or 'NA'}"),
        "rescued_reads_total": totals,
    }


# ------------------------------------------------------ input alignments (X13 / E3)


def _qc_rows(params: dict) -> list[dict] | None:
    root = params.get("qc_summary_dir")
    if not root or root in {"null", "None"}:
        return None
    path = Path(root) / "cmp.qc.summary.tsv"
    if not path.is_file():
        return None
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def input_alignments_summary(params: dict) -> dict | None:
    """Cohort view of the BAM/CRAM input provenance (scripts/header_provenance.py columns of the
    QC summary): how many inputs per reference fingerprint, per aligner configuration, per @PG
    chain and per upstream unmapped filter. Counts only -- no sample ID, path or read group."""
    rows = [r for r in (_qc_rows(params) or [])
            if r.get("input_ref_fingerprint") not in (None, "", "NA")]
    if not rows:
        return None

    def tally(key):
        out: dict[str, int] = defaultdict(int)
        for r in rows:
            out[key(r)] += 1
        return dict(sorted(out.items(), key=lambda kv: (-kv[1], kv[0])))

    def value(r, column):
        return r.get(column) or "NA"

    def reference(r):
        return (f"{value(r, 'input_ref_fingerprint')[:12]} {value(r, 'input_ref_build')}, "
                f"{value(r, 'input_ref_contigs')} contigs (ALT {value(r, 'input_ref_alt_contigs')}, "
                f"decoy {value(r, 'input_ref_decoy_contigs')}, HLA {value(r, 'input_ref_hla_contigs')}, "
                f"viral {value(r, 'input_ref_viral_contigs')}, chrEBV {value(r, 'input_ref_has_ebv')})")

    def aligner(r):
        name = f"{value(r, 'upstream_aligner')} {value(r, 'upstream_aligner_version')}"
        options = r.get("input_bwa_options")
        return f"{name} {options}" if options not in (None, "", "NA") else name

    return {
        "inputs": len(rows),
        "references": tally(reference),
        "aligners": tally(aligner),
        "pg_chains": tally(lambda r: value(r, "input_pg_chain")),
        "unmapped_filter": tally(lambda r: value(r, "input_unmapped_filter")),
    }


def render_input_alignments(summary: dict | None) -> list[str]:
    if not summary:
        return []
    refs, aligners = summary["references"], summary["aligners"]
    out = ["INPUT ALIGNMENTS (BAM/CRAM headers; per sample: QC_SUMMARY input_* columns)",
           f"  {summary['inputs']} input(s): {len(refs)} reference fingerprint(s); "
           f"{len(aligners)} aligner configuration(s)"]
    out.append("  reference (md5 of @SQ SN:LN, first 12 hex)")
    out += [f"    {n:>6}  {key}" for key, n in refs.items()]
    out.append("  aligner (from @PG; bwa mem options reduced to letters and numbers)")
    out += [f"    {n:>6}  {key}" for key, n in aligners.items()]
    out.append("  @PG chain (program:version)")
    out += [f"    {n:>6}  {key}" for key, n in summary["pg_chains"].items()]
    out.append("  upstream unmapped filter")
    out += [f"    {n:>6}  {key}" for key, n in summary["unmapped_filter"].items()]
    if len(refs) > 1:
        out.append("  *** more than one input reference: unmapped reads were defined against different")
        out.append("      contig sets (decoy / viral / ALT), so extraction is not equivalent across inputs.")
    out.append("")
    return out


# ------------------------------------------------------------ database versions (X14)


def _sibling(name: str):
    """A helper module from scripts/, loaded by path (this file is also loaded by path in tests)."""
    import importlib.util
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"cmp_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DATABASE_PARAMS = ("hg38_db", "t2t_phix_db", "pangenome_db", "cram_reference", "kraken_db", "metaphlan_db",
                   "humann3_nucleotide_db", "humann3_protein_db", "validation_ref", "validation_taxdb")


def database_versions(params: dict, stages: list[str] | None, manifest: Path | None = None) -> dict:
    """conf/databases.tsv compared with the databases this run used (scripts/database_versions.py).
    Never raises: a provenance failure must not mask the run's own outcome."""
    try:
        dv = _sibling("database_versions")
        paths = {name: params.get(name) for name in DATABASE_PARAMS}
        paths = {k: (None if v in (None, "", "null", "None") else v) for k, v in paths.items()}
        # The taxonomy DB alignment validation reads when no prebuilt reference is given (main.nf).
        if not paths.get("validation_ref") and not paths.get("validation_taxdb") and paths.get("kraken_db"):
            paths["validation_taxdb"] = f"{paths['kraken_db']}/taxDB"
        summary = dv.record(paths, stages, manifest or dv.MANIFEST)
        summary["refseq_snapshot"] = params.get("validation_refseq_snapshot") or None
        return summary
    except Exception as exc:   # noqa: BLE001 -- recorded, never fatal
        return {"error": f"{type(exc).__name__}: {exc}", "rows": []}


def render_database_versions(summary: dict | None) -> list[str]:
    if not summary:
        return []
    if summary.get("error"):
        return ["DATABASE VERSIONS", f"  could not be recorded: {summary['error']}", ""]
    lines = _sibling("database_versions").render(summary)
    if summary.get("refseq_snapshot") and lines:
        lines.insert(-1, f"  {'validation RefSeq':<22}snapshot {summary['refseq_snapshot']} (--validation_refseq_snapshot)")
    return lines


# ------------------------------------------------------- sample outcomes (X05/X06)


def sample_outcomes(params: dict) -> dict | None:
    """How the samples ended, from the QC summary's terminal_state column: counts per state and
    the gated-stage statuses. Counts only -- no sample ID enters the record (per-sample states stay
    in cmp.qc.summary.tsv on the cluster). None when there is no QC summary or no such column."""
    root = params.get("qc_summary_dir")
    if not root or root in {"null", "None"}:
        return None
    path = Path(root) / "cmp.qc.summary.tsv"
    if not path.is_file():
        return None
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "terminal_state" not in rows[0]:
        return None

    def tally(column):
        out: dict[str, int] = defaultdict(int)
        for row in rows:
            out[row.get(column) or "NA"] += 1
        return dict(sorted(out.items()))

    by_state = tally("terminal_state")
    return {
        "samples": len(rows),
        "by_state": by_state,
        "failed": sum(n for state, n in by_state.items() if state.startswith("failed_at_")),
        # QC floor (review #34): not microbial negatives; excluded from the cohort tables unless
        # --qc_floor_action flag.
        "not_assessable": sum(n for state, n in by_state.items() if state.startswith("not_assessable_")),
        "qc_floor_action": params.get("qc_floor_action") or "exclude",
        "metaphlan_status": tally("metaphlan_status") if "metaphlan_status" in rows[0] else {},
        "humann_status": tally("humann_status") if "humann_status" in rows[0] else {},
    }


def render_sample_outcomes(outcomes: dict | None) -> list[str]:
    if not outcomes:
        return []
    out = ["SAMPLE OUTCOMES (QC_SUMMARY/cmp.qc.summary.tsv, terminal_state)"]
    for state, n in outcomes["by_state"].items():
        out.append(f"  {n:>6}  {state}")
    if outcomes.get("not_assessable"):
        action = outcomes.get("qc_floor_action", "exclude")
        out.append(f"  *** {outcomes['not_assessable']} sample(s) NOT ASSESSABLE (QC floor): no unmapped reads to")
        out.append("      profile, unmapped reads filtered upstream, or too few reads after host depletion.")
        out.append("      They are NOT microbial negatives. " + (
            "--qc_floor_action exclude: they are in no cohort table." if action == "exclude" else
            "--qc_floor_action flag: they WERE profiled and are in the cohort tables; drop them before analysis."))
    if outcomes["failed"]:
        out.append(f"  *** {outcomes['failed']} sample(s) failed at a per-sample stage: they are missing")
        out.append("      from the cohort tables because a task failed, not because they hold no microbes.")
        out.append("      Which samples and stages: the terminal_state column of the QC summary.")
    for column in ("metaphlan_status", "humann_status"):
        if outcomes.get(column):
            out.append(f"  {column:<18}" + "   ".join(f"{k} {v}" for k, v in outcomes[column].items()))
    out.append("")
    return out


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

    # Which release each database was, and whether it is still the frozen one (X14).
    for line in render_database_versions(record.get("databases")):
        add(line)

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

    # What each task actually ran (K06/K07): the CONDA ENVS section above is what was specified.
    for line in render_tool_versions((record.get("environment") or {}).get("tool_versions") or {}):
        add(line)

    rescue = (record.get("outputs") or {}).get("low_score_rescue")
    if rescue:
        add("LOW-SCORE RESCUE (BAM/CRAM, --extract_low_score_mapped)")
        add(f"  inputs        {rescue['inputs']}   " +
            "   ".join(f"{mode} {n}" for mode, n in rescue["by_mode"].items()))
        for reason, n in rescue["by_reason"].items():
            add(f"  {n:>6}  {reason}")
        for aligner, n in rescue["by_aligner"].items():
            add(f"  {n:>6}  upstream {aligner}")
        totals = rescue["rescued_reads_total"]
        for column in LOW_AS_COUNTS:
            value = totals.get(column)
            add(f"  {column:<30}{'-' if value is None else value:>14}  ({totals.get(column + '_samples', 0)} samples)")
        add("")

    for line in render_input_alignments((record.get("outputs") or {}).get("input_alignments")):
        add(line)

    for line in render_sample_outcomes((record.get("outputs") or {}).get("sample_outcomes")):
        add(line)

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
    # strict | warn; default: the launch record's --tool_version_check, else strict.
    parser.add_argument("--tool-version-check", choices=("strict", "warn"), default=None)
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
        # --extract_low_score_mapped, aggregated over the BAM/CRAM inputs (counts only; L8).
        "low_score_rescue": low_score_rescue_summary(launch_params),
        # BAM/CRAM header provenance: reference fingerprints, aligner chains (counts only; X13).
        "input_alignments": input_alignments_summary(launch_params),
        # How each sample ended (terminal_state), counts only (audit X05/X06).
        "sample_outcomes": sample_outcomes(launch_params),
    }

    # Database releases and fingerprints, compared with conf/databases.tsv (X14).
    record["databases"] = database_versions(
        launch_params, (record.get("invocation") or {}).get("stages_requested"))

    environment = record.setdefault("environment", {})
    envs = conda_envs(environment.get("conda_envs") or {}, launch_params.get("conqur_sha"))
    environment["conda_envs_used"] = envs
    software = software_summary(envs)
    if software:
        environment["software"] = software

    # Runtime tool versions, read back from each task's work dir (K06/K07).
    check = args.tool_version_check or launch_params.get("tool_version_check") or "strict"
    if check in ("null", "None"):
        check = "strict"
    environment["tool_versions"] = collect_tool_versions(
        rows, (record.get("invocation") or {}).get("work_dir"),
        Path(__file__).resolve().parents[1], check)

    target.write_text(json.dumps(record, indent=2) + "\n")

    report = render(record, processes, software, args.tracedir, args.session)
    report_path = args.tracedir / "RUN_REPORT.txt"
    with report_path.open("a") as handle:
        handle.write(report)
    print(f"Run report: {report_path}")

    tool_problems = environment["tool_versions"]["problems"]
    for entry in tool_problems:
        print(f"[provenance] TOOL VERSION MISMATCH {entry['status']}: {entry['detail']} "
              f"({', '.join(entry['processes'])}; {entry['tasks']} task(s))")
    # strict: a non-zero exit, which the onComplete handler reports as an error. MISSING and
    # INSTALLED_MISMATCH can only be found here, after the tasks ran.
    blocking = [e for e in tool_problems if not e.get("waived")]
    return 1 if blocking and check == "strict" else 0


if __name__ == "__main__":
    raise SystemExit(main())
