#!/usr/bin/env python3
"""Database versions: what each reference database a run read actually was (audit X14).

"Path is not a version." lib/Provenance.groovy's data() digest (relative path, size and mtime to
depth 3, truncated to 12 hex) is a -resume cache key, not something a reader can compare with a
release: KrakenUniq `krakenUniq_8_8_2023`, MetaPhlAn `mpa_vJun23`, ChocoPhlAn/UniRef and the host
indexes could each be rebuilt or replaced in place and every run before and after would print the
same line. conf/databases.tsv freezes each database's release label and a cheap fingerprint; every
run recomputes the fingerprint for the databases it used and records whether it matches (finalize
-> RUN_REPORT "DATABASE VERSIONS", JSON `databases`). Ported from pksProfiler c08098d
(conf/databases.tsv, scripts/record_database_versions.sh); the listing fingerprint is computed the
same way, so the two pipelines' frozen values for the shared files are directly comparable.

FINGERPRINTS, and why they are cheap. KrakenUniq's database.kdb alone is 574 GB; hashing it on every
run is tens of minutes of Lustre reads to learn nothing. So, per database:

  sha256:listing   sha256 over the sorted "<relative path>\\t<size in bytes>\\n" lines of every
                   regular file (symlinks not followed) -- or "<name>\\t<size>\\n" for one large file.
                   O(file count); changes when a database is rebuilt, re-downloaded or extended.
  sha256:content   sha256 of the file itself, for a single file at or below CMP_DB_HASH_MAX_BYTES
                   (default 64 MiB).
  key (key_files)  a real content check of small release-defining files, recomputed every run:
                   `sha256:files` = sha256 over sorted "<name>\\t<sha256 of content>\\n" lines of the
                   files matching the key globs (each at most CMP_DB_KEY_MAX_BYTES, default 256 MiB);
                   `mmi:dictionary` = sha256 of a minimap2 index's header + sequence dictionary (the
                   first few kB: build settings and every sequence name and length;
                   pksProfiler scripts/describe_minimap2_index.py).
  newest_mtime     the newest file modification time, recorded and compared but NOT part of the
                   identity: an rsync or restore rewrites mtimes without changing the data.

What none of this catches: an in-place byte edit that keeps every size and every key file. The
integrity check is `sha256sum` of the big files, run deliberately, not per run.

Verdicts: frozen | frozen (path differs) | CHANGED: listing | CHANGED: key | not_in_manifest |
absent | unused (stage not run) | unset. A record, not a gate: it never fails a run.

    database_versions.py record --manifest conf/databases.tsv --stages a,b NAME=PATH ... [--json OUT]
    database_versions.py freeze NAME=PATH ...    # print manifest rows to paste, for a new freeze

Standard library only. Reads file metadata, a few small key files and the first kB of .mmi files.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import json
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path

MANIFEST = Path(__file__).resolve().parents[1] / "conf" / "databases.tsv"
MANIFEST_FIELDS = ["name", "param", "stage", "database", "release", "path", "select", "size_bytes",
                   "file_count", "fingerprint_method", "fingerprint", "newest_mtime", "key_method",
                   "key_files", "key_sha256", "frozen_on", "source"]
RECORD_FIELDS = ["name", "param", "path", "kind", "status", "release", "size_bytes", "file_count",
                 "fingerprint_method", "fingerprint", "newest_mtime", "key_method", "key_sha256",
                 "verdict", "note"]

CONTENT_MAX = int(os.environ.get("CMP_DB_HASH_MAX_BYTES", 64 * 1024 ** 2))
KEY_MAX = int(os.environ.get("CMP_DB_KEY_MAX_BYTES", 256 * 1024 ** 2))
MAX_FILES = 500_000

# Each database parameter, the stage that reads it (main.nf / RunRecord.stages names), and the
# rows a run records for it. bracken_kmer_distrib is a view of kraken_db: Bracken has no database
# of its own; its read-length models are the database<N>mers.kmer_distrib files inside it.
PARAM_STAGE = {
    "hg38_db": "host_depletion", "t2t_phix_db": "host_depletion", "pangenome_db": "pangenome_depletion",
    "cram_reference": "host_depletion",
    "kraken_db": "classification", "metaphlan_db": "classification",
    "humann3_nucleotide_db": "humann3", "humann3_protein_db": "humann3",
    "validation_ref": "alignment_validation", "validation_taxdb": "alignment_validation",
}
DERIVED = {"bracken_kmer_distrib": ("kraken_db", "database*mers.kmer_distrib")}


# --------------------------------------------------------------------------- fingerprints


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _iso(mtime: float | None) -> str:
    if mtime is None:
        return ""
    return datetime.fromtimestamp(int(mtime), timezone.utc).isoformat()


def _regular_files(root: Path, select: str):
    """(relative path, size, mtime) of every regular file below root matching `select`."""
    out = []
    for dirpath, dirnames, filenames in os.walk(root):          # followlinks=False
        dirnames.sort()
        for name in filenames:
            full = Path(dirpath) / name
            if full.is_symlink() or not full.is_file():
                continue
            rel = full.relative_to(root).as_posix()
            if select and select != "*" and not fnmatch.fnmatchcase(rel, select):
                continue
            st = full.stat()
            out.append((rel, st.st_size, st.st_mtime))
            if len(out) > MAX_FILES:
                raise RuntimeError(f"more than {MAX_FILES} files")
    return sorted(out)


def listing_sha256(entries) -> str:
    text = "".join(f"{rel}\t{size}\n" for rel, size, *_ in sorted(entries))
    return hashlib.sha256(text.encode("utf-8", "surrogateescape")).hexdigest()


def mmi_dictionary_sha256(path: Path) -> str:
    """sha256 of a minimap2 .mmi's header + sequence dictionary (pksProfiler describe_minimap2_index)."""
    with path.open("rb") as handle:
        if handle.read(4) != b"MMI\x02":
            raise ValueError(f"{path.name}: not a minimap2 index")
        _w, _k, _b, n_seq, _flag = struct.unpack("<5I", handle.read(20))
        for _ in range(n_seq):
            raw = handle.read(1)
            if not raw:
                raise ValueError(f"{path.name}: truncated sequence dictionary")
            handle.read(raw[0])
            if len(handle.read(4)) != 4:
                raise ValueError(f"{path.name}: truncated sequence dictionary")
        size = handle.tell()
    with path.open("rb") as handle:
        return hashlib.sha256(handle.read(size)).hexdigest()


def key_fingerprint(path: Path, kind: str, method: str, globs: str, files) -> str:
    """The key check (see the module docstring). '' when the row declares none."""
    if not method:
        return ""
    if method == "mmi:dictionary":
        target = path if kind == "file" else None
        if target is None:
            raise ValueError("mmi:dictionary needs a single .mmi file")
        return mmi_dictionary_sha256(target)
    if method != "sha256:files":
        raise ValueError(f"unknown key method {method!r}")
    patterns = [g.strip() for g in (globs or "").split(",") if g.strip()]
    if kind == "file":
        chosen = [(path.name, path)]
    else:
        chosen = [(rel, path / rel) for rel, *_ in files
                  if any(fnmatch.fnmatchcase(rel, g) for g in patterns)]
    if not chosen:
        raise ValueError(f"no key file matches {globs!r}")
    lines = []
    for rel, full in sorted(chosen):
        if full.stat().st_size > KEY_MAX:
            raise ValueError(f"key file {rel} is above CMP_DB_KEY_MAX_BYTES")
        lines.append(f"{rel}\t{sha256_file(full)}\n")
    return hashlib.sha256("".join(lines).encode()).hexdigest()


def fingerprint(name: str, param: str, path_text: str | None, select: str = "*",
                key_method: str = "", key_files: str = "") -> dict:
    """One RECORD_FIELDS row (without verdict) for the database at `path_text`."""
    row = {f: "" for f in RECORD_FIELDS} | {"name": name, "param": param, "path": path_text or ""}
    if not path_text or path_text in ("null", "None"):
        row.update(kind="unset", status="unset")
        return row
    path = Path(path_text)
    try:
        if path.is_dir():
            files = _regular_files(path, select)
            row["kind"] = "directory"
            if not files:
                row.update(status="unreadable", note="no regular files")
                return row
            row.update(fingerprint_method="sha256:listing", fingerprint=listing_sha256(files),
                       file_count=len(files), size_bytes=sum(s for _, s, _ in files),
                       newest_mtime=_iso(max(m for *_, m in files)))
        elif path.is_file():
            st = path.stat()
            files = [(path.name, st.st_size, st.st_mtime)]
            row.update(kind="file", file_count=1, size_bytes=st.st_size, newest_mtime=_iso(st.st_mtime))
            if st.st_size <= CONTENT_MAX:
                row.update(fingerprint_method="sha256:content", fingerprint=sha256_file(path))
            else:
                row.update(fingerprint_method="sha256:listing", fingerprint=listing_sha256(files))
        else:
            # A Bowtie2 / minimap2 index given as a prefix: the siblings sharing it.
            parent = path.parent
            files = []
            if parent.is_dir():
                files = sorted((p.name, p.stat().st_size, p.stat().st_mtime) for p in parent.iterdir()
                               if p.name.startswith(path.name) and p.is_file() and not p.is_symlink())
            if not files:
                row.update(kind="missing", status="absent")
                return row
            row.update(kind="prefix", fingerprint_method="sha256:listing", fingerprint=listing_sha256(files),
                       file_count=len(files), size_bytes=sum(s for _, s, _ in files),
                       newest_mtime=_iso(max(m for *_, m in files)))
            path = parent
        row["status"] = "ok"
        if key_method:
            row["key_method"] = key_method
            row["key_sha256"] = key_fingerprint(path, row["kind"], key_method, key_files, files)
    except (OSError, ValueError, RuntimeError) as exc:
        row.update(status="unreadable", note=str(exc)[:120])
    return row


# --------------------------------------------------------------------------- manifest


def read_manifest(path: Path = MANIFEST) -> dict[str, dict]:
    rows = {}
    if not Path(path).is_file():
        return rows
    with Path(path).open() as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    for record in csv.DictReader(lines, delimiter="\t"):
        if record.get("name"):
            if record["name"] in rows:
                raise ValueError(f"{path}: duplicate database {record['name']}")
            rows[record["name"]] = record
    return rows


def verdict(row: dict, frozen: dict | None) -> tuple[str, str]:
    if row["status"] in ("unset", "unused"):
        return row["status"], row.get("note", "")
    if row["status"] == "absent":
        return "absent", ""
    if row["status"] != "ok":
        return "unreadable", row.get("note", "")
    if not frozen:
        return "not_in_manifest", ""
    if row["fingerprint"] != frozen.get("fingerprint"):
        return "CHANGED: listing", f"frozen {frozen.get('fingerprint', '')[:16]} on {frozen.get('frozen_on', '?')}"
    if frozen.get("key_sha256") and row.get("key_sha256") != frozen["key_sha256"]:
        return "CHANGED: key", f"frozen key {frozen['key_sha256'][:16]} on {frozen.get('frozen_on', '?')}"
    notes = []
    if frozen.get("newest_mtime") and row.get("newest_mtime") != frozen["newest_mtime"]:
        notes.append(f"newest mtime {row.get('newest_mtime')} differs from frozen {frozen['newest_mtime']}"
                     " (copied or touched; sizes and keys match)")
    if frozen.get("path") and os.path.normpath(frozen["path"]) != os.path.normpath(row["path"]):
        return "frozen (path differs)", "; ".join(notes)
    return "frozen", "; ".join(notes)


def record(paths: dict[str, str | None], stages: list[str] | None, manifest: Path = MANIFEST) -> dict:
    """Fingerprint every database parameter in `paths` (param -> path) for the stages this run
    requested, compare with the manifest, and return {manifest, manifest_sha256, rows}."""
    frozen = read_manifest(manifest)
    stages = set(stages) if stages is not None else None
    wanted = [(p, p, "*") for p in PARAM_STAGE if p in paths] + \
             [(n, p, sel) for n, (p, sel) in DERIVED.items() if p in paths]
    rows = []
    for name, param, select in wanted:
        stage = PARAM_STAGE[param]
        entry = frozen.get(name, {})
        if stages is not None and stage not in stages:
            row = {f: "" for f in RECORD_FIELDS} | {"name": name, "param": param,
                                                    "path": paths.get(param) or "", "status": "unused",
                                                    "note": f"stage {stage} not requested"}
        else:
            row = fingerprint(name, param, paths.get(param), entry.get("select") or select,
                              entry.get("key_method", ""), entry.get("key_files", ""))
        row["release"] = entry.get("release", "")
        row["verdict"], note = verdict(row, entry or None)
        row["note"] = "; ".join(dict.fromkeys(n for n in (row.get("note"), note) if n))
        rows.append(row)
    return {
        "manifest": str(manifest),
        "manifest_sha256": sha256_file(Path(manifest)) if Path(manifest).is_file() else None,
        "rows": rows,
    }


def human(n) -> str:
    try:
        value = float(n)
    except (TypeError, ValueError):
        return "-"
    for unit, size in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("KB", 1e3)):
        if value >= size:
            return f"{value / size:.1f} {unit}"
    return f"{value:.0f} B"


def render(summary: dict | None) -> list[str]:
    """RUN_REPORT's DATABASE VERSIONS section."""
    if not summary:
        return []
    out = ["DATABASE VERSIONS (conf/databases.tsv; listing fingerprint + key files, recomputed this run)"]
    rows = [r for r in summary.get("rows", []) if r.get("status") not in ("unset",)]
    if not rows:
        out += ["  no database parameter set", ""]
        return out
    for r in rows:
        if r["status"] == "unused":
            out.append(f"  {r['name']:<22}unused ({r.get('note', '')})")
            continue
        out.append(f"  {r['name']:<22}{r['verdict']}")
        out.append(f"  {'':<22}release  {r.get('release') or '(not in manifest)'}")
        if r.get("fingerprint"):
            out.append(f"  {'':<22}{r['fingerprint_method']} {r['fingerprint'][:16]}  "
                       f"{human(r.get('size_bytes'))} in {r.get('file_count')} file(s)"
                       + (f"; key {r['key_method']} {r['key_sha256'][:16]}" if r.get("key_sha256") else ""))
        if r.get("note"):
            out.append(f"  {'':<22}{r['note']}")
    changed = [r["name"] for r in rows if str(r.get("verdict", "")).startswith("CHANGED")]
    if changed:
        out.append(f"  *** {', '.join(changed)}: NOT the frozen database. Results are not comparable with runs")
        out.append("      against the frozen release until the change is understood and conf/databases.tsv updated.")
    out.append("")
    return out


# --------------------------------------------------------------------------- CLI


def _pairs(items: list[str]) -> dict[str, str]:
    out = {}
    for item in items:
        name, sep, value = item.partition("=")
        if not sep:
            raise SystemExit(f"expected NAME=PATH, got {item!r}")
        out[name] = value
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("record")
    r.add_argument("--manifest", type=Path, default=MANIFEST)
    r.add_argument("--stages", default=None, help="comma-separated stages requested (default: all)")
    r.add_argument("--json", type=Path, default=None)
    r.add_argument("pairs", nargs="*")
    f = sub.add_parser("freeze", help="print manifest rows (fingerprint columns) for NAME=PATH")
    f.add_argument("--manifest", type=Path, default=MANIFEST)
    f.add_argument("pairs", nargs="*")
    a = ap.parse_args()

    paths = _pairs(a.pairs)
    if a.cmd == "record":
        stages = [s for s in a.stages.split(",") if s] if a.stages is not None else None
        summary = record(paths, stages, a.manifest)
        if a.json:
            a.json.write_text(json.dumps(summary, indent=2) + "\n")
        print("\n".join(render(summary)))
        return 0

    frozen = read_manifest(a.manifest)
    today = datetime.now(timezone.utc).date().isoformat()
    writer = csv.DictWriter(sys.stdout, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
    writer.writeheader()
    for name, path in paths.items():
        entry = dict(frozen.get(name, {}))
        param, select = DERIVED.get(name, (name, "*"))
        row = fingerprint(name, param, path, entry.get("select") or select,
                          entry.get("key_method", ""), entry.get("key_files", ""))
        entry.update({k: row[k] for k in ("size_bytes", "file_count", "fingerprint_method", "fingerprint",
                                          "newest_mtime", "key_sha256")})
        entry.update(name=name, param=param, path=path, frozen_on=today,
                     stage=entry.get("stage") or PARAM_STAGE.get(param, ""), select=entry.get("select") or select)
        writer.writerow(entry)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
