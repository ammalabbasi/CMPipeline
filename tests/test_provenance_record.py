#!/usr/bin/env python3
"""Checks for scripts/finalize_provenance.py.

The completion half of the run record is driven from workflow.onComplete, which is
awkward to exercise end to end. The parsing and rendering it delegates to are plain
functions, so they are tested here directly: a failed run must say `failed`, a dirty
working tree must be called out, and a header-only trace must degrade gracefully
rather than crash the handler.

Also (audit C06/D08): the record is one file per launch (<sessionId>_<runName>.json),
it carries SHA-256 digests of the cohort tables and nothing of their content, and it
lists only the conda envs this run used -- lockfile digests, and for a pre-built prefix
env its conda-meta and the ConQuR/cqrReg DESCRIPTIONs -- never a walk of the cache.
All inputs here are synthetic.
"""

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "finalize_provenance", ROOT / "scripts" / "finalize_provenance.py"
)
FP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FP)

TRACE_HEADER = (
    "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t"
    "%cpu\tpeak_rss\tpeak_vmem\trchar\twchar"
)


def trace_row(name, status="COMPLETED", peak="1.3 GB"):
    return (
        f"1\tab/cdef12\t1000\t{name}\t{status}\t0\t2026-09-22 10:00:00\t"
        f"5m\t4m\t180.0%\t{peak}\t2.0 GB\t1 GB\t1 GB"
    )


class ByteParsingTests(unittest.TestCase):
    def test_parses_units(self):
        self.assertAlmostEqual(FP.to_bytes("1.3 GB"), 1.3e9)
        self.assertAlmostEqual(FP.to_bytes("512 MB"), 512e6)

    def test_missing_values_are_none(self):
        for text in ("-", "", "   "):
            self.assertIsNone(FP.to_bytes(text))

    def test_human_bytes_roundtrip(self):
        self.assertEqual(FP.human_bytes(1.3e9), "1.30 GB")
        self.assertEqual(FP.human_bytes(None), "-")


class TraceParsingTests(unittest.TestCase):
    def test_header_only_trace_yields_no_rows(self):
        """The three trace-*.txt files already in the repo are exactly this shape."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "trace.txt"
            path.write_text(TRACE_HEADER + "\n")
            self.assertEqual(FP.read_trace(path), [])

    def test_missing_trace_yields_no_rows(self):
        self.assertEqual(FP.read_trace(Path("/nonexistent/trace.txt")), [])

    def test_per_sample_tag_is_stripped(self):
        self.assertEqual(FP.process_name("filterReads (SAMPLE1)"), "filterReads")
        self.assertEqual(FP.process_name("consensus_taxa"), "consensus_taxa")

    def test_summary_groups_by_process_and_counts_failures(self):
        rows = FP.read_trace_rows = [
            dict(zip(TRACE_HEADER.split("\t"), trace_row("filterReads (A)").split("\t"))),
            dict(zip(TRACE_HEADER.split("\t"), trace_row("filterReads (B)", peak="2.6 GB").split("\t"))),
            dict(zip(TRACE_HEADER.split("\t"), trace_row("mapReads (A)", status="FAILED").split("\t"))),
        ]
        summary = {entry["process"]: entry for entry in FP.summarise_processes(rows)}
        self.assertEqual(summary["filterReads"]["tasks"], 2)
        self.assertEqual(summary["filterReads"]["failed"], 0)
        self.assertAlmostEqual(summary["filterReads"]["peak_rss_max"], 2.6e9)
        self.assertEqual(summary["mapReads"]["failed"], 1)


class RenderingTests(unittest.TestCase):
    def base_record(self, dirty=False, status="completed"):
        return {
            "record_version": 1,
            "status": status,
            "code": {
                "commit": "abc123",
                "commit_source": "git",
                "dirty": dirty,
                "modified_files": ["M main.nf"] if dirty else [],
                "scripts_digest": "aaa/content",
                "modules_digest": "bbb/content",
            },
            "invocation": {"stages_requested": ["host_depletion"], "work_dir": "/w"},
            "dependencies": {"hg38_db": "ccc/meta"},
            "environment": {"nextflow_version": "26.04.6"},
            "run": {"started": "t0", "task_counts": {"succeeded": "3"}},
        }

    def test_dirty_tree_is_called_out(self):
        text = FP.render(self.base_record(dirty=True), [], {}, Path("/t"), "s1")
        self.assertIn("UNCOMMITTED CHANGES", text)
        self.assertIn("NOT attributable", text)

    def test_clean_tree_is_not_flagged(self):
        text = FP.render(self.base_record(dirty=False), [], {}, Path("/t"), "s1")
        self.assertNotIn("UNCOMMITTED CHANGES", text)
        self.assertIn("working tree  clean", text)

    def test_failed_run_says_failed(self):
        record = self.base_record(status="failed")
        record["run"]["error_message"] = "Process mapReads terminated"
        text = FP.render(record, [], {}, Path("/t"), "s1")
        self.assertIn("FAILED", text)
        self.assertIn("Process mapReads terminated", text)

    def test_empty_trace_explains_itself(self):
        text = FP.render(self.base_record(), [], {}, Path("/t"), "s1")
        self.assertIn("per-process detail unavailable", text)

    def test_process_table_is_rendered_when_present(self):
        processes = [
            {"process": "filterReads", "tasks": 2, "failed": 0,
             "peak_rss_max": 1.3e9, "peak_rss_median": 1.2e9,
             "realtime_sample": "5m", "pct_cpu_sample": "180%"},
        ]
        text = FP.render(self.base_record(), processes, {}, Path("/t"), "s1")
        self.assertIn("filterReads", text)
        self.assertIn("1.30 GB", text)
        self.assertIn("sized from", text)


class FinaliseTests(unittest.TestCase):
    def test_missing_launch_record_is_not_fatal(self):
        """A provenance failure must never mask the run's own outcome."""
        with tempfile.TemporaryDirectory() as tmp:
            import sys
            argv = sys.argv
            sys.argv = ["finalize_provenance.py", "--tracedir", tmp, "--session", "nope"]
            try:
                self.assertEqual(FP.main(), 0)
            finally:
                sys.argv = argv

    def test_status_flips_to_failed(self):
        with tempfile.TemporaryDirectory() as tmp:
            tracedir = Path(tmp)
            (tracedir / "provenance").mkdir()
            record = {"status": "started", "code": {}, "invocation": {},
                      "dependencies": {}, "environment": {}, "run": {}}
            (tracedir / "provenance" / "sess-1_happy_turing.json").write_text(json.dumps(record))
            import sys
            argv = sys.argv
            sys.argv = ["finalize_provenance.py", "--tracedir", str(tracedir),
                        "--session", "sess-1_happy_turing", "--success", "false", "--exit-status", "1"]
            try:
                FP.main()
            finally:
                sys.argv = argv
            written = json.loads((tracedir / "provenance" / "sess-1_happy_turing.json").read_text())
            self.assertEqual(written["status"], "failed")
            self.assertTrue((tracedir / "RUN_REPORT.txt").is_file())

    def test_two_launches_of_one_session_keep_two_records(self):
        """-resume reuses the session id; the run name keeps each launch's record (C06)."""
        with tempfile.TemporaryDirectory() as tmp:
            tracedir = Path(tmp)
            (tracedir / "provenance").mkdir()
            for run_name, flag in (("first_run", "false"), ("second_run", "true")):
                record = {"status": "started", "invocation": {"params": {"run_decontam": flag}},
                          "environment": {}, "run": {}}
                (tracedir / "provenance" / f"sess-1_{run_name}.json").write_text(json.dumps(record))
                run_main(["--tracedir", str(tracedir), "--session", f"sess-1_{run_name}",
                          "--success", "true"])
            first = json.loads((tracedir / "provenance" / "sess-1_first_run.json").read_text())
            second = json.loads((tracedir / "provenance" / "sess-1_second_run.json").read_text())
            self.assertEqual(first["invocation"]["params"]["run_decontam"], "false")
            self.assertEqual(second["invocation"]["params"]["run_decontam"], "true")


def run_main(args):
    import sys
    argv = sys.argv
    sys.argv = ["finalize_provenance.py", *args]
    try:
        return FP.main()
    finally:
        sys.argv = argv


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# A marker that must never appear in the record or the report: proves only digests leave.
MARKER = "SYNTHETIC_CONTENT_MARKER_7f3a"


def make_outdir(root: Path) -> dict:
    """A synthetic results tree: cohort tables plus per-sample files that must be skipped."""
    out = root / "RESULTS"
    files = {
        "CONSENSUS_TAXA/bracken.metaphlan.common.genus.mpa.report.txt": b"t\tA\n" + MARKER.encode(),
        "CONSENSUS_TAXA/consensus_status.tsv": b"x\n",
        "BRACKEN/bracken.genus.mpa.report.txt": b"g\n",
        "04_DECONTAMINATION/genus/genus.decontam_pkg_decontaminated.csv": b"a,b\n",
        "04_DECONTAMINATION/genus/genus.contaminants_final.tsv": b"c\n",
        "05_BATCH_CORRECTION/genus/corrected/conqur_corrected.tsv": b"d\n",
        "05_BATCH_CORRECTION/genus/batch_correction_sample_status.tsv": b"e\n",
        "06_VALIDATION/validation_mask_report.tsv": b"f\n",
        "QC_SUMMARY/cmp.qc.summary.tsv": b"q\n",
        # must NOT be hashed: per-sample outputs and read files
        "BRACKEN/S1.bracken.G.report.txt": b"per-sample\n",
        "06_VALIDATION/per_sample/S1.validation_stats.tsv": b"per-sample\n",
        "06_VALIDATION/S1.validated_reads.fasta.gz": b"reads\n",
        "CONSENSUS_TAXA/stray.fastq.gz": b"reads\n",
    }
    for rel, data in files.items():
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        (out / rel).write_bytes(data)
    params = {
        "outdir": str(out),
        "consensus_taxa_dir": str(out / "CONSENSUS_TAXA"),
        "krakenuniq_bracken_dir": str(out / "BRACKEN"),
        "metaphlan4_dir": str(out / "METAPHLAN4"),       # absent: skipped quietly
        "humann3_dir": None,
        "decontam_dir": str(out / "04_DECONTAMINATION"),
        "batch_corr_dir": str(out / "05_BATCH_CORRECTION"),
        "validation_dir": str(out / "06_VALIDATION"),
        "qc_summary_dir": str(out / "QC_SUMMARY"),
    }
    return {"params": params, "files": files}


class OutputChecksumTests(unittest.TestCase):
    def test_cohort_tables_are_hashed_and_per_sample_files_are_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree = make_outdir(Path(tmp))
            result = FP.output_checksums(tree["params"])
            by_path = {entry["path"]: entry for entry in result["files"]}
            expected = {rel for rel in tree["files"]
                        if "S1" not in rel and "fastq" not in rel}
            self.assertEqual(set(by_path), expected)
            for rel in expected:
                self.assertEqual(by_path[rel]["sha256"], sha(tree["files"][rel]))
                self.assertEqual(by_path[rel]["bytes"], len(tree["files"][rel]))
            self.assertFalse(result["truncated"])

    def test_a_changed_table_changes_its_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree = make_outdir(Path(tmp))
            before = FP.output_checksums(tree["params"])
            target = Path(tree["params"]["qc_summary_dir"]) / "cmp.qc.summary.tsv"
            target.write_bytes(b"q2\n")
            after = FP.output_checksums(tree["params"])
            pick = lambda r: {e["path"]: e["sha256"] for e in r["files"]}["QC_SUMMARY/cmp.qc.summary.tsv"]
            self.assertNotEqual(pick(before), pick(after))

    def test_no_outdir_yields_empty_list(self):
        self.assertEqual(FP.output_checksums({"outdir": "/nonexistent/RESULTS"})["files"], [])

    def test_record_and_report_carry_digests_not_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree = make_outdir(Path(tmp))
            tracedir = Path(tmp) / "pipeline_info"
            (tracedir / "provenance").mkdir(parents=True)
            record = {"status": "started", "invocation": {"params": tree["params"]},
                      "environment": {}, "run": {}}
            (tracedir / "provenance" / "s_r.json").write_text(json.dumps(record))
            run_main(["--tracedir", str(tracedir), "--session", "s_r", "--success", "true"])
            text = (tracedir / "provenance" / "s_r.json").read_text()
            report = (tracedir / "RUN_REPORT.txt").read_text()
            written = json.loads(text)
            self.assertEqual(written["outputs"]["cohort_tables"]["file_count"], 9)
            self.assertNotIn(MARKER, text)
            self.assertNotIn(MARKER, report)
            self.assertIn("COHORT TABLES", report)


LOCK_LINES = [
    "https://conda.anaconda.org/conda-forge/linux-64/samtools-1.21-h96c455f_1.conda#aaaa",
    "https://conda.anaconda.org/conda-forge/linux-64/htslib-1.21-h566b1c6_1.tar.bz2#bbbb",
    "https://conda.anaconda.org/conda-forge/linux-64/libzlib-1.3.1-hb9d3cd8_2.conda#cccc",
]


def make_lock(repo: Path, name="samtools_env") -> Path:
    lock = repo / "conda_envs" / "lock" / f"{name}.linux-64.txt"
    lock.parent.mkdir(parents=True, exist_ok=True)
    lock.write_text("# platform: linux-64\n@EXPLICIT\n" + "\n".join(LOCK_LINES) + "\n")
    (repo / "conda_envs" / f"{name}.yml").write_text("name: x\ndependencies:\n  - samtools=1.21\n")
    return lock


def make_prefix(root: Path, lines=LOCK_LINES, conqur_sha="ff233085") -> Path:
    """A fake built env: conda-meta for each lock line, an R library with a conda-owned
    package and two installed after the solve (ConQuR from GitHub, cqrReg from CRAN)."""
    prefix = root / "batch_env"
    meta = prefix / "conda-meta"
    meta.mkdir(parents=True)
    for line in lines:
        url, md5 = line.split("#")
        fname = url.rsplit("/", 1)[-1]
        m = FP.PACKAGE_FILE.match(fname)
        (meta / fname.replace(".conda", ".json").replace(".tar.bz2", ".json")).write_text(json.dumps({
            "name": m["name"], "version": m["version"], "build": m["build"],
            "url": url, "md5": md5, "files": [],
        }))
    (meta / "r-vegan-2.6_4-r43h0_0.json").write_text(json.dumps({
        "name": "r-vegan", "version": "2.6_4", "build": "r43h0_0",
        "files": ["lib/R/library/vegan/DESCRIPTION", "lib/R/library/vegan/R/vegan"],
    }))
    lib = prefix / "lib" / "R" / "library"
    for pkg, desc in {
        "vegan": "Package: vegan\nVersion: 2.6-4\n",
        "ConQuR": f"Package: ConQuR\nVersion: 2.0\nRemoteType: github\nRemoteUsername: ivartb\n"
                  f"RemoteRepo: ConQuR_par\nRemoteSha: {conqur_sha}\nDescription: long\n  continued\n",
        "cqrReg": "Package: cqrReg\nVersion: 1.2.1\nRepository: CRAN\n",
    }.items():
        (lib / pkg).mkdir(parents=True)
        (lib / pkg / "DESCRIPTION").write_text(desc)
    return prefix


class CondaEnvTests(unittest.TestCase):
    def test_lockfile_env_is_digested_and_mapped_to_its_yml(self):
        with tempfile.TemporaryDirectory() as tmp:
            lock = make_lock(Path(tmp))
            entry = FP.env_record(str(lock))
            self.assertEqual(entry["kind"], "lockfile")
            self.assertEqual(entry["sha256"], sha(lock.read_bytes()))
            self.assertTrue(entry["source_yml"].endswith("conda_envs/samtools_env.yml"))
            self.assertIn("samtools=1.21=h96c455f_1", entry["packages"])
            self.assertIn("htslib=1.21=h566b1c6_1", entry["packages"])

    def test_one_changed_lock_line_changes_the_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            lock = make_lock(Path(tmp))
            before = FP.env_record(str(lock))
            lock.write_text(lock.read_text().replace("samtools-1.21", "samtools-1.22"))
            after = FP.env_record(str(lock))
            self.assertNotEqual(before["sha256"], after["sha256"])
            self.assertNotEqual(before["explicit_sha256"], after["explicit_sha256"])

    def test_prefix_env_reads_conda_meta_and_r_descriptions(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = make_prefix(Path(tmp))
            entry = FP.env_record(str(prefix), conqur_sha="ff233085")
            self.assertEqual(entry["kind"], "prefix")
            self.assertEqual(entry["r_packages"]["ConQuR"]["Version"], "2.0")
            self.assertEqual(entry["r_packages"]["ConQuR"]["RemoteSha"], "ff233085")
            self.assertEqual(entry["r_packages"]["cqrReg"]["Version"], "1.2.1")
            self.assertNotIn("Description", entry["r_packages"]["ConQuR"])
            self.assertEqual(entry["r_packages_not_from_conda"],
                             {"ConQuR": "2.0@ff233085", "cqrReg": "1.2.1"})
            self.assertTrue(entry["conqur_sha_matches_param"])
            self.assertIn("r-vegan=2.6_4=r43h0_0", entry["packages"])

    def test_conqur_mismatch_is_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = make_prefix(Path(tmp), conqur_sha="0000000")
            entry = FP.env_record(str(prefix), conqur_sha="ff233085")
            self.assertFalse(entry["conqur_sha_matches_param"])

    def test_env_built_from_a_lock_has_the_locks_explicit_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            lock = make_lock(Path(tmp))
            prefix = make_prefix(Path(tmp))
            self.assertEqual(FP.env_record(str(lock))["explicit_sha256"],
                             FP.env_record(str(prefix))["explicit_sha256"])

    def test_only_this_runs_envs_are_listed_and_the_cache_is_not_walked(self):
        """D08 acceptance: a fake cache with a current and a stale env; only the
        recorded spec is reported, named by its parameter, with ConQuR 2.0 ff233085."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lock = make_lock(root)
            prefix = make_prefix(root)
            cache = root / "conda_cache"
            for stale in ("env-stale000", "env-current1"):
                (cache / stale / "conda-meta").mkdir(parents=True)
                (cache / stale / "conda-meta" / "samtools-1.9-h0_0.json").write_text("{}")
            tracedir = root / "pipeline_info"
            (tracedir / "provenance").mkdir(parents=True)
            record = {"status": "started", "run": {},
                      "invocation": {"params": {"conqur_sha": "ff233085"}},
                      "environment": {"conda_envs": {"samtools_env": str(lock),
                                                     "batch_corr_env": str(prefix)}}}
            (tracedir / "provenance" / "s_r.json").write_text(json.dumps(record))
            run_main(["--tracedir", str(tracedir), "--session", "s_r", "--success", "true",
                      "--conda-cache-dir", str(cache)])
            written = json.loads((tracedir / "provenance" / "s_r.json").read_text())
            envs = written["environment"]["conda_envs_used"]
            self.assertEqual(set(envs), {"samtools_env", "batch_corr_env"})
            text = json.dumps(written)
            self.assertNotIn("env-stale000", text)
            self.assertNotIn("samtools=1.9", text)
            self.assertIn("ConQuR=2.0@ff233085", written["environment"]["software"]["batch_corr_env"])
            report = (tracedir / "RUN_REPORT.txt").read_text()
            self.assertIn("CONDA ENVS THIS RUN USED", report)
            self.assertIn("2 R package(s) not installed by conda", report)

    def test_missing_spec_is_recorded_as_absent(self):
        self.assertEqual(FP.env_record("/nonexistent/lock.linux-64.txt")["kind"], "absent")


class LaunchHalfTests(unittest.TestCase):
    """The launch half lives in Groovy (lib/, main.nf). These pin its source so the C06/D08
    fields cannot silently disappear; the -preview check in TEST_NOW.sh exercises it."""

    MAIN = (ROOT / "main.nf").read_text()
    RUNRECORD = (ROOT / "lib" / "RunRecord.groovy").read_text()
    PROVENANCE = (ROOT / "lib" / "Provenance.groovy").read_text()

    def test_record_key_is_session_plus_run_name(self):
        self.assertIn('RunRecord.path(params.tracedir, "${workflow.sessionId}_${workflow.runName}")', self.MAIN)
        self.assertIn("'--session',     \"${wfMeta.sessionId}_${wfMeta.runName}\"", self.MAIN)

    def test_sheet_and_metadata_are_content_hashed(self):
        self.assertIn("digest: Provenance.content(params.sample)", self.MAIN)
        self.assertIn("digest: Provenance.content(params.metadata_file)", self.MAIN)
        self.assertIn('return "sha256:${digest.digest().encodeHex()}"', self.PROVENANCE)

    def test_validation_dependencies_and_envs_are_recorded(self):
        for key in ("validation_ref", "validation_taxdb", "validation_groups",
                    "validation_cache_dir", "envs_digest", "conda_envs"):
            self.assertRegex(self.MAIN, rf"\b{key}\s*:", key)
        self.assertNotIn("'--conda-cache-dir', condaCache", self.MAIN)

    def test_new_stages_are_listed(self):
        self.assertIn("'alignment_validation'", self.RUNRECORD)
        self.assertIn("'cohort_report'", self.RUNRECORD)


if __name__ == "__main__":
    unittest.main()
