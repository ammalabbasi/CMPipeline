#!/usr/bin/env python3
"""Runtime tool versions and the task-environment guard (audit K06/K07/X04; ported from
pksProfiler F16: 26e196d, 2bc7ca9, 2c0f8c8).

CRC50 attempt 1 (Job 12380606) ran minimap2 2.24 from ~/.bashrc's tools/minimap2 in every mapReads
task while the lockfile pinned 2.28, and RUN_REPORT -- which lists what the lockfiles specify --
said nothing. These tests cover the three parts of the fix:

  1. scripts/record_tool_versions.sh, executed with FAKE tools on PATH: ok, a version off its pin,
     a tool shadowed by an earlier PATH entry (the CRC50 case; must fail strict mode), a tool the
     env lacks, a shell function, no env, retries, and pins read from a lockfile, a yml and a
     prefix's conda-meta;
  2. static wiring: every process that runs a tool from a conda env records it, each recorded
     tool's package is in that env's lockfile, the beforeScript guard, exit 78 -> terminate, the
     trace's workdir column, the param and the onComplete handler;
  3. scripts/finalize_provenance.py reading the records back (TOOL VERSIONS / TOOL VERSION
     MISMATCH, MISSING, INSTALLED_MISMATCH, strict exit status);

and, when nextflow and conda are available, 4. a Nextflow run of a one-process pipeline against a
fake conda env with a fake "login" minimap2 earlier on PATH: the repo's beforeScript makes the env
win; without it the task stops with exit 78 and the run terminates.

Every input is synthetic: fake tools print a version string and nothing else.
"""

import csv
import importlib.util
import json
import os
import re
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "record_tool_versions.sh"
MODULES = sorted((REPO / "Modules").glob("*.nf"))
LOCK = REPO / "conda_envs" / "lock"

SPEC = importlib.util.spec_from_file_location("finalize_provenance", REPO / "scripts" / "finalize_provenance.py")
FP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FP)

# Processes with a conda env that deliberately do not record. antiSMASH belongs to the separate
# run_antismash.nf entry point (its own params; not part of main.nf).
EXEMPT = {"antiSMASH": "separate entry point run_antismash.nf, not part of main.nf"}

GUARD = ('unset LD_LIBRARY_PATH LIBRARY_PATH CPATH PYTHONHOME PERL5LIB R_LIBS R_LIBS_USER R_LIBS_SITE; '
         'if command -v conda >/dev/null 2>&1; then export PATH="$(conda info --base)/bin:$PATH"; fi')

LOCK_HEADER = "# platform: linux-64\n@EXPLICIT\n"


def lock_line(name, version, build="h0_0", ext="conda", channel="bioconda"):
    return f"https://conda.anaconda.org/{channel}/linux-64/{name}-{version}-{build}.{ext}#0123456789abcdef0123456789abcdef\n"


def fake_tool(directory, name, output, exit_code=0, stderr=False):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    redirect = " >&2" if stderr else ""
    path.write_text(f"#!/bin/sh\necho '{output}'{redirect}\nexit {exit_code}\n")
    path.chmod(0o755)
    return path


def clean_env(**extra):
    """A minimal environment: system tools only, plus what the test adds."""
    env = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp"), "LANG": "C",
           "CMP_TOOL_VERSION_BACKOFF": "0"}
    env.update(extra)
    return env


def run_record(mode, env_spec, tools, env, process="mapReads", out_dir=None, script=SCRIPT):
    out = Path(out_dir) / ".tool_versions.tsv"
    result = subprocess.run(["bash", str(script), mode, str(env_spec), process, str(out), *tools],
                            env=env, capture_output=True, text=True)
    rows = []
    if out.is_file():
        with out.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    return result, {r["tool"]: r for r in rows}


class RecordScript(unittest.TestCase):
    """scripts/record_tool_versions.sh with fake tools on PATH."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="cmp_toolver_"))
        self.prefix = self.tmp / "env"
        (self.prefix / "conda-meta").mkdir(parents=True)
        self.lock = self.tmp / "minimap2_env.linux-64.txt"
        self.lock.write_text(LOCK_HEADER + lock_line("minimap2", "2.28", "h577a1d6_4", "tar.bz2")
                             + lock_line("samtools", "1.21") + lock_line("htslib", "1.24")
                             + lock_line("r-base", "4.3.3", "hc28aa72_21", channel="conda-forge"))
        fake_tool(self.prefix / "bin", "minimap2", "2.28-r1209")
        fake_tool(self.prefix / "bin", "samtools", "samtools 1.21")
        fake_tool(self.prefix / "bin", "bgzip", "bgzip (htslib) 1.24")
        self.login = self.tmp / "login_tools"

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def env(self, path_head=(), **extra):
        path = ":".join([*map(str, path_head), str(self.prefix / "bin"), "/usr/bin", "/bin"])
        return clean_env(PATH=path, CONDA_PREFIX=str(self.prefix), **extra)

    def test_pinned_tools_from_the_env_are_ok(self):
        result, rows = run_record("strict", self.lock, ["minimap2", "samtools", "bgzip:htslib"],
                                  self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({t: r["status"] for t, r in rows.items()},
                         {"minimap2": "ok", "samtools": "ok", "bgzip": "ok"})
        self.assertEqual(rows["minimap2"]["version"], "2.28")
        self.assertEqual(rows["minimap2"]["pinned"], "2.28")      # synthetic lockfile (self.lock), not the repo's
        self.assertEqual(rows["bgzip"]["package"], "htslib")
        self.assertEqual(rows["bgzip"]["pinned"], "1.24")
        self.assertEqual(rows["minimap2"]["path"], str(self.prefix / "bin" / "minimap2"))
        self.assertEqual(rows["minimap2"]["attempts"], "1")

    def test_record_columns(self):
        _result, _rows = run_record("strict", self.lock, ["minimap2"], self.env(), out_dir=self.tmp)
        header = (self.tmp / ".tool_versions.tsv").read_text().splitlines()[0].split("\t")
        self.assertEqual(header, ["process", "tool", "package", "status", "version", "pinned", "path",
                                  "conda_prefix", "env", "raw", "attempts", "read_exit", "check"])

    def test_login_minimap2_ahead_of_the_env_fails_strict(self):
        """The CRC50 Job 12380606 case: tools/minimap2 2.24 earlier on PATH than the env."""
        fake_tool(self.login, "minimap2", "2.24-r1122")
        result, rows = run_record("strict", self.lock, ["minimap2", "samtools"],
                                  self.env(path_head=[self.login]), out_dir=self.tmp)
        self.assertEqual(result.returncode, 78)
        self.assertEqual(rows["minimap2"]["status"], "OUTSIDE_ENV")
        self.assertEqual(rows["minimap2"]["version"], "2.24")
        self.assertEqual(rows["minimap2"]["path"], str(self.login / "minimap2"))
        self.assertEqual(rows["samtools"]["status"], "ok")
        self.assertIn("outside this task's conda environment", result.stderr)
        self.assertIn("pinned minimap2=2.28", result.stderr)

    def test_shadowing_is_recorded_but_not_fatal_in_warn_mode(self):
        fake_tool(self.login, "minimap2", "2.24-r1122")
        result, rows = run_record("warn", self.lock, ["minimap2"],
                                  self.env(path_head=[self.login]), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(rows["minimap2"]["status"], "OUTSIDE_ENV")
        self.assertIn("WARNING", result.stderr)

    def test_shadowing_with_the_same_version_still_fails(self):
        """Path, not just version: an identical version from another install is still not the env."""
        fake_tool(self.login, "samtools", "samtools 1.21")
        result, rows = run_record("strict", self.lock, ["samtools"],
                                  self.env(path_head=[self.login]), out_dir=self.tmp)
        self.assertEqual(result.returncode, 78)
        self.assertEqual(rows["samtools"]["status"], "OUTSIDE_ENV")

    def test_tool_the_env_lacks_resolves_to_the_login_path(self):
        fake_tool(self.login, "fastp", "fastp 0.20.1")
        result, rows = run_record("strict", self.lock, ["fastp"], self.env(path_head=[]) | {
            "PATH": f"{self.prefix / 'bin'}:{self.login}:/usr/bin:/bin"}, out_dir=self.tmp)
        self.assertEqual(result.returncode, 78)
        self.assertEqual(rows["fastp"]["status"], "OUTSIDE_ENV")

    def test_version_off_its_pin(self):
        fake_tool(self.prefix / "bin", "minimap2", "2.24-r1122")
        result, rows = run_record("strict", self.lock, ["minimap2"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 78)
        self.assertEqual(rows["minimap2"]["status"], "VERSION_MISMATCH")
        result, rows = run_record("warn", self.lock, ["minimap2"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(rows["minimap2"]["status"], "VERSION_MISMATCH")

    def test_more_specific_release_satisfies_the_pin(self):
        lock = self.tmp / "short.txt"
        lock.write_text(LOCK_HEADER + lock_line("minimap2", "2.28"))
        fake_tool(self.prefix / "bin", "minimap2", "2.28.1")
        result, rows = run_record("strict", lock, ["minimap2"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows["minimap2"]["status"], "ok")
        # ...but 2.281 is not 2.28.
        fake_tool(self.prefix / "bin", "minimap2", "2.281")
        result, rows = run_record("strict", lock, ["minimap2"], self.env(), out_dir=self.tmp)
        self.assertEqual(rows["minimap2"]["status"], "VERSION_MISMATCH")

    def test_not_found(self):
        result, rows = run_record("strict", self.lock, ["krakenuniq"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 78)
        self.assertEqual(rows["krakenuniq"]["status"], "NOT_FOUND")

    def test_no_env_active(self):
        env = clean_env(PATH=f"{self.prefix / 'bin'}:/usr/bin:/bin")
        result, rows = run_record("strict", self.lock, ["minimap2"], env, out_dir=self.tmp)
        self.assertEqual(result.returncode, 78)
        self.assertEqual(rows["minimap2"]["status"], "NO_ENV")

    def test_exported_shell_function_shadows_the_env(self):
        """A function outranks every PATH entry; the integration test's krakenuniq stub is one."""
        result = subprocess.run(
            ["bash", "-c", 'minimap2() { echo "2.28-r1209"; }; export -f minimap2; '
                           f'bash "{SCRIPT}" strict "{self.lock}" mapReads "{self.tmp}/.tool_versions.tsv" minimap2'],
            env=self.env(), capture_output=True, text=True)
        self.assertEqual(result.returncode, 78)
        with (self.tmp / ".tool_versions.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["status"], "OUTSIDE_ENV")
        self.assertEqual(rows[0]["path"], "(function minimap2)")

    def test_unpinned_package(self):
        fake_tool(self.prefix / "bin", "pigz", "pigz 2.8")
        result, rows = run_record("strict", self.lock, ["pigz"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(rows["pigz"]["status"], "unpinned")

    def test_package_named_after_the_colon(self):
        fake_tool(self.prefix / "bin", "Rscript", "Rscript (R) version 4.3.3 (2024-02-29)")
        result, rows = run_record("strict", self.lock, ["Rscript:r-base"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((rows["Rscript"]["package"], rows["Rscript"]["pinned"], rows["Rscript"]["status"]),
                         ("r-base", "4.3.3", "ok"))

    def test_hyphenated_package_names_do_not_collide(self):
        """r-base-4.3.3 must not be read as package r with version base."""
        lock = self.tmp / "r.txt"
        lock.write_text(LOCK_HEADER + lock_line("r-base", "4.3.3") + lock_line("python-dateutil", "2.9.0")
                        + lock_line("python", "3.11.16"))
        fake_tool(self.prefix / "bin", "python", "Python 3.11.16")
        result, rows = run_record("strict", lock, ["python"], self.env(), out_dir=self.tmp)
        self.assertEqual(rows["python"]["pinned"], "3.11.16")
        self.assertEqual(rows["python"]["status"], "ok")

    def test_yml_pins(self):
        yml = self.tmp / "env.yml"
        yml.write_text("dependencies:\n  - bioconda::minimap2=2.28  # pinned\n  - samtools>=1.0\n")
        fake_tool(self.prefix / "bin", "minimap2", "2.24-r1122")
        result, rows = run_record("strict", yml, ["minimap2", "samtools"], self.env(), out_dir=self.tmp)
        self.assertEqual(rows["minimap2"]["status"], "VERSION_MISMATCH")
        self.assertEqual(rows["samtools"]["status"], "unpinned")   # a range is not a pin
        self.assertEqual(result.returncode, 78)

    def test_prefix_pins_from_conda_meta(self):
        """--batch_corr_env <prefix>: the prefix's own conda-meta gives the versions."""
        (self.prefix / "conda-meta" / "r-base-4.3.3-hc28aa72_21.json").write_text("{}")
        fake_tool(self.prefix / "bin", "Rscript", "Rscript (R) version 4.3.3 (2024-02-29)")
        result, rows = run_record("strict", self.prefix, ["Rscript:r-base"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows["Rscript"]["pinned"], "4.3.3")

    def test_version_on_stderr_and_paths_stripped(self):
        fake_tool(self.prefix / "bin", "samtools", "samtools 1.21", stderr=True)
        # bowtie2 prints its own path, which holds digits that must not read as a version.
        fake_tool(self.prefix / "bin", "bowtie2",
                  f"{self.prefix}/bin/python3.10/bowtie2-align-s version 2.5.4")
        lock = self.tmp / "b.txt"
        lock.write_text(LOCK_HEADER + lock_line("bowtie2", "2.5.4") + lock_line("samtools", "1.21"))
        result, rows = run_record("strict", lock, ["samtools", "bowtie2"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows["bowtie2"]["version"], "2.5.4")
        self.assertEqual(rows["samtools"]["version"], "1.21")

    def test_a_failed_read_is_retried(self):
        state = self.tmp / "calls"
        tool = self.prefix / "bin" / "minimap2"
        tool.write_text(f"#!/bin/sh\necho x >> '{state}'\n"
                        f"[ $(wc -l < '{state}') -ge 2 ] && echo 2.28-r1209 && exit 0\nexit 1\n")
        tool.chmod(0o755)
        result, rows = run_record("strict", self.lock, ["minimap2"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((rows["minimap2"]["status"], rows["minimap2"]["attempts"]), ("ok", "2"))

    def test_unreadable_version_is_recorded_not_fatal(self):
        fake_tool(self.prefix / "bin", "samtools", "no version here")
        result, rows = run_record("strict", self.lock, ["samtools"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(rows["samtools"]["status"], "version_unread")
        self.assertEqual(rows["samtools"]["attempts"], "1")   # exit 0 with output: not retried
        fake_tool(self.prefix / "bin", "samtools", "", exit_code=3)
        result, rows = run_record("strict", self.lock, ["samtools"],
                                  self.env(CMP_TOOL_VERSION_ATTEMPTS="2"), out_dir=self.tmp)
        self.assertEqual(result.returncode, 0)
        self.assertEqual((rows["samtools"]["attempts"], rows["samtools"]["read_exit"]), ("2", "3"))
        self.assertIn("exit 3 after 2 attempt(s)", rows["samtools"]["raw"])

    def test_alias_is_accepted_for_exactly_that_build(self):
        scripts = self.tmp / "scripts"
        scripts.mkdir()
        shutil.copy(SCRIPT, scripts / SCRIPT.name)
        (scripts / "tool_version_aliases.tsv").write_text("bracken\t3.1p1\t3.0.1\twhy\n")
        lock = self.tmp / "k.txt"
        lock.write_text(LOCK_HEADER + lock_line("bracken", "3.1p1"))
        fake_tool(self.prefix / "bin", "bracken", "Bracken v3.0.1")
        result, rows = run_record("strict", lock, ["bracken"], self.env(), out_dir=self.tmp,
                                  script=scripts / SCRIPT.name)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows["bracken"]["status"], "ok")

    def test_bad_mode(self):
        result, _rows = run_record("loose", self.lock, ["minimap2"], self.env(), out_dir=self.tmp)
        self.assertEqual(result.returncode, 2)

    def test_real_lockfiles_parse(self):
        """Every lockfile yields exact pins for the tools the pipeline records from it."""
        result, rows = run_record("warn", LOCK / "minimap2_env.linux-64.txt", ["minimap2", "samtools"],
                                  self.env(), out_dir=self.tmp)
        # Expected pin comes from the env spec, so a tool update (2.28 -> 2.31, 2026-10-26) needs no test edit.
        spec = (REPO / "conda_envs" / "minimap2_env.yml").read_text()
        self.assertEqual(rows["minimap2"]["pinned"], re.search(r"minimap2=([\d.]+)", spec).group(1))
        self.assertEqual(rows["samtools"]["pinned"], FP.env_pins(str(LOCK / "minimap2_env.linux-64.txt"))["samtools"])


# ------------------------------------------------------------------------------ static wiring

PROCESS_RE = re.compile(r"^process\s+(\w+)\s*\{", re.M)
RECORD_RE = re.compile(r"ToolVersions\.record\(task, params, '([^']+)'\)")
CONDA_RE = re.compile(r'^\s*conda\s+"\$\{params\.(\w+)\}"', re.M)
PARAM_ENV_RE = re.compile(r"params\.(\w+_env)\s*=\s*\"\$\{projectDir\}/([^\"]+)\"")


def processes():
    """{process name: (module file, body text)}."""
    found = {}
    for module in MODULES:
        text = module.read_text()
        starts = list(PROCESS_RE.finditer(text))
        for index, match in enumerate(starts):
            end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
            found[match.group(1)] = (module, text[match.start():end])
    return found


def env_lockfiles():
    main = (REPO / "main.nf").read_text()
    envs = {name: REPO / rel for name, rel in PARAM_ENV_RE.findall(main)}
    envs["batch_corr_env"] = LOCK / "batch_correction_env.linux-64.txt"   # BatchEnv.lockfile's source
    return envs


class Wiring(unittest.TestCase):

    def test_every_conda_process_records_its_tools(self):
        missing = []
        for name, (module, body) in processes().items():
            if CONDA_RE.search(body) and not RECORD_RE.search(body) and name not in EXEMPT:
                missing.append(f"{module.name}:{name}")
        self.assertEqual(missing, [], "processes running a conda env without ToolVersions.record")

    def test_exemptions_are_not_stale(self):
        for name in EXEMPT:
            self.assertIn(name, processes())
            self.assertNotIn("ToolVersions.record", processes()[name][1])

    def test_every_process_has_a_conda_env(self):
        """The record's OUTSIDE_ENV check assumes every task activates its own env."""
        bare = [n for n, (_m, body) in processes().items() if not CONDA_RE.search(body)]
        self.assertEqual(bare, [])

    def test_recorded_packages_are_in_the_process_env(self):
        envs = env_lockfiles()
        problems = []
        for name, (module, body) in processes().items():
            records = RECORD_RE.findall(body)
            if not records:
                continue
            self.assertEqual(len(records), 1, f"{name}: one record line per process")
            env_param = CONDA_RE.search(body).group(1)
            pins = FP.env_pins(str(envs[env_param]))
            self.assertTrue(pins, f"{name}: no pins parsed from {envs[env_param]}")
            for entry in records[0].split():
                tool, _, package = entry.partition(":")
                package = package or tool
                if package not in pins:
                    problems.append(f"{name}: {tool} ({package}) not in {envs[env_param].name}")
        self.assertEqual(problems, [])

    def test_python_tools_name_their_package(self):
        """`python3` in a record maps to the python package in both the script and the finaliser."""
        for name, (_module, body) in processes().items():
            for records in RECORD_RE.findall(body):
                for entry in records.split():
                    tool, _, package = entry.partition(":")
                    if tool == "python3":
                        # Without :python, record_tool_versions.sh looks for a package called
                        # python3, finds no pin and records `unpinned`: no version check at all.
                        self.assertEqual(package, "python", f"{name}: record python3 as python3:python")

    def test_record_line_runs_before_the_tools(self):
        """The record is the first command of each script it is in (bar shell settings and vars)."""
        for name, (_module, body) in processes().items():
            if "${ToolVersions.record" not in body:
                continue
            # From the opening quote of the script string that holds the record line.
            before = body.split("${ToolVersions.record", 1)[0].rsplit('"""', 1)[1]
            commands = [l.strip() for l in before.splitlines()
                        if l.strip() and not l.strip().startswith("#")]
            for line in commands:
                self.assertRegex(line, r"^(set -|export |[A-Z_0-9]+=|out[0-9]=)",
                                 f"{name}: '{line}' runs before the tool-version record")

    def test_python_tool_matches_the_interpreter_the_script_calls(self):
        for name, (_module, body) in processes().items():
            records = RECORD_RE.findall(body)
            if not records:
                continue
            script = body.split('"""', 1)[1]
            tools = [entry.split(":")[0] for entry in records[0].split()]
            calls = set(re.findall(r"^\s*(python3?)\s", script, re.M))
            for call in calls:
                self.assertIn(call, tools, f"{name} calls {call} but records {tools}")

    def test_before_script_guard(self):
        config = (REPO / "nextflow.config").read_text()
        match = re.search(r"^process\.beforeScript\s*=\s*'(.*)'\s*$", config, re.M)
        self.assertIsNotNone(match)
        self.assertEqual(match.group(1), GUARD)
        integration = (REPO / "tests" / "integration" / "integration.config").read_text()
        self.assertIn(f"beforeScript = '{GUARD}' +", integration,
                      "integration.config's Bracken beforeScript must start with the same guard")
        # PYTHONNOUSERSITE / PYTHONPATH stay in the env scope (exported after activation).
        self.assertRegex(config, r"PYTHONNOUSERSITE\s*=\s*'1'")

    def test_exit_78_terminates_in_every_error_strategy(self):
        base = (REPO / "conf" / "base.config").read_text()
        strategies = re.findall(r"errorStrategy\s*=\s*\{([^}]*)\}", base)
        self.assertGreaterEqual(len(strategies), 2)
        for strategy in strategies:
            self.assertTrue(strategy.strip().startswith("task.exitStatus == 78 ? 'terminate'"), strategy)
        self.assertIn("exit 78", SCRIPT.read_text())
        self.assertIn("EXIT_STATUS = 78", (REPO / "lib" / "ToolVersions.groovy").read_text())

    def test_no_module_overrides_the_strategy_for_recorded_processes(self):
        """A process-level errorStrategy would win over base.config's 78 -> terminate."""
        for name, (_module, body) in processes().items():
            if RECORD_RE.search(body) and "errorStrategy" in body:
                # Only FASTQC carries one, and its per_sample label's config selector overrides it.
                self.assertIn("label 'per_sample'", body, name)

    def test_trace_has_workdir(self):
        config = (REPO / "nextflow.config").read_text()
        fields = re.search(r"trace\.fields\s*=\s*'([^']+)'", config).group(1).split(",")
        for needed in ("name", "status", "hash", "workdir", "attempt", "cpus"):
            self.assertIn(needed, fields)

    def test_param_and_handler(self):
        main = (REPO / "main.nf").read_text()
        self.assertIn("params.tool_version_check = 'strict'", main)
        self.assertIn("ToolVersions.validate(params.tool_version_check)", main)
        self.assertIn("'env', '-u', 'LD_LIBRARY_PATH'", main)
        self.assertIn("proc.exitValue() != 0", main)

    def test_finaliser_knows_every_recording_process(self):
        names = FP.recording_processes(REPO)
        for alias in ("FASTQC1", "FASTQCHG38", "mapReads", "Bracken", "BatchCorrection", "qcSummary"):
            self.assertIn(alias, names)
        self.assertNotIn("antiSMASH", names)

    def test_documented(self):
        params_doc = (REPO / "docs" / "parameters.md").read_text()
        self.assertIn("--tool_version_check", params_doc)
        self.assertIn("TOOL VERSIONS", (REPO / "docs" / "outputs.md").read_text())
        self.assertIn("tool_version_check", (REPO / "docs" / "CHANGELOG.md").read_text())
        self.assertIn("tool_version_check", (REPO / "docs" / "hpc.md").read_text())


# -------------------------------------------------------------------------------- finaliser

TRACE_HEADER = ["task_id", "hash", "native_id", "name", "status", "exit", "workdir"]
RECORD_HEADER = ["process", "tool", "package", "status", "version", "pinned", "path", "conda_prefix",
                 "env", "raw", "attempts", "read_exit", "check"]


class Finaliser(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="cmp_toolver_fp_"))
        self.lock = self.tmp / "minimap2_env.linux-64.txt"
        self.lock.write_text(LOCK_HEADER + lock_line("minimap2", "2.28") + lock_line("samtools", "1.21"))
        self.prefix = self.tmp / "env"
        (self.prefix / "conda-meta").mkdir(parents=True)
        for name in ("minimap2-2.28-h0_0.json", "samtools-1.21-h0_0.json"):
            (self.prefix / "conda-meta" / name).write_text("{}")
        self.rows = []

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def task(self, name, records, status="COMPLETED"):
        workdir = self.tmp / "work" / f"{len(self.rows):02d}" / "abcdef"
        workdir.mkdir(parents=True)
        if records is not None:
            with (workdir / ".tool_versions.tsv").open("w") as handle:
                handle.write("\t".join(RECORD_HEADER) + "\n")
                for record in records:
                    handle.write("\t".join(record.get(c, "") for c in RECORD_HEADER) + "\n")
        self.rows.append({"name": name, "status": status, "workdir": str(workdir), "hash": "ab/cdef12"})

    def record(self, tool, version, status="ok", path=None, prefix=None, package=None, env=None):
        prefix = str(prefix or self.prefix)
        return {"process": "mapReads", "tool": tool, "package": package or tool, "status": status,
                "version": version, "pinned": "", "path": path or f"{prefix}/bin/{tool}",
                "conda_prefix": prefix, "env": str(env or self.lock), "attempts": "1", "read_exit": "0"}

    def collect(self, mode="strict"):
        return FP.collect_tool_versions(self.rows, None, REPO, mode)

    def test_clean_run(self):
        self.task("mapReads (S1)", [self.record("minimap2", "2.28"), self.record("samtools", "1.21")])
        self.task("mapReads (S2)", [self.record("minimap2", "2.28"), self.record("samtools", "1.21")])
        summary = self.collect()
        self.assertEqual(summary["problems"], [])
        self.assertEqual(summary["tasks_recorded"], 2)
        minimap2 = [t for t in summary["tools"] if t["tool"] == "minimap2"][0]
        self.assertEqual((minimap2["tasks"], minimap2["pinned"], minimap2["installed"]), (2, "2.28", "2.28"))
        text = "\n".join(FP.render_tool_versions(summary))
        self.assertIn("TOOL VERSIONS (as run)", text)
        self.assertNotIn("TOOL VERSION MISMATCH", text)

    def test_warn_run_with_a_shadowed_tool_is_reported(self):
        self.task("mapReads (S1)", [self.record("minimap2", "2.24", status="OUTSIDE_ENV",
                                                path="/home/u/tools/minimap2/minimap2")])
        summary = self.collect("warn")
        self.assertEqual([p["status"] for p in summary["problems"]], ["OUTSIDE_ENV"])
        text = "\n".join(FP.render_tool_versions(summary))
        self.assertIn("TOOL VERSION MISMATCH", text)
        self.assertIn("were NOT produced by the pinned tools", text)
        self.assertIn("/home/u/tools/minimap2/minimap2", text)

    def test_task_level_warn_is_listed_but_not_blocking(self):
        """2026-10-26: the integration test relaxes one task (ext.tool_version_check = 'warn') for its
        krakenuniq stub; a strict run must list that row but not fail on it alone."""
        waived = dict(self.record("krakenuniq", "1.0.4", status="OUTSIDE_ENV", path="(function krakenuniq)"),
                      check="warn")
        self.task("Bracken (S1)", [waived])
        summary = self.collect("strict")
        self.assertEqual([(p["status"], p["waived"]) for p in summary["problems"]], [("OUTSIDE_ENV", True)])
        self.assertIn("listed, not blocking", "\n".join(FP.render_tool_versions(summary)))
        # the same problem from a strict task is blocking
        self.task("Bracken (S2)", [dict(waived, check="strict")])
        self.assertEqual([p["waived"] for p in self.collect("strict")["problems"]], [False])

    def test_bracken_alias_is_accepted(self):
        """bioconda bracken 3.1p1 prints v3.0.1 (scripts/tool_version_aliases.tsv)."""
        lock = self.tmp / "k.txt"
        lock.write_text(LOCK_HEADER + lock_line("bracken", "3.1p1"))
        self.task("Bracken (S1)", [self.record("bracken", "3.0.1", env=lock)])
        self.assertEqual(self.collect("strict")["problems"], [])

    def test_spec_is_rechecked(self):
        """A record the task called ok is re-judged against its spec's pins here."""
        self.task("mapReads (S1)", [self.record("minimap2", "2.24")])
        statuses = [p["status"] for p in self.collect()["problems"]]
        self.assertIn("VERSION_MISMATCH", statuses)

    def test_installed_mismatch(self):
        (self.prefix / "conda-meta" / "samtools-1.21-h0_0.json").rename(
            self.prefix / "conda-meta" / "samtools-1.19-h0_0.json")
        lock = self.tmp / "unpinned.yml"
        lock.write_text("dependencies:\n  - samtools\n")
        self.task("mapReads (S1)", [self.record("samtools", "1.21", status="unpinned", env=lock)])
        self.assertEqual([p["status"] for p in self.collect()["problems"]], ["INSTALLED_MISMATCH"])

    def test_missing_record(self):
        self.task("mapReads (S1)", None)
        self.task("FASTQC1 (S1)", None)
        self.task("mapReads (S2)", None, status="FAILED")   # a failed task need not have one
        problems = self.collect()["problems"]
        self.assertEqual(sorted(p["processes"][0] for p in problems), ["FASTQC1", "mapReads"])
        self.assertTrue(all(p["status"] == "MISSING" for p in problems))

    def test_unreadable_work_dir_is_counted(self):
        self.rows.append({"name": "mapReads (S9)", "status": "CACHED", "workdir": str(self.tmp / "gone")})
        summary = self.collect()
        self.assertEqual(summary["tasks_unreadable"], 1)
        self.assertIn("could not be checked", "\n".join(FP.render_tool_versions(summary)))

    def test_no_sample_ids_in_the_summary(self):
        self.task("mapReads (SAMPLE_X1)", [self.record("minimap2", "2.24", status="OUTSIDE_ENV",
                                                       path="/x/minimap2")])
        self.assertNotIn("SAMPLE_X1", json.dumps(self.collect("warn")))

    def test_main_exit_status_and_report(self):
        tracedir = self.tmp / "pipeline_info"
        (tracedir / "provenance").mkdir(parents=True)
        self.task("mapReads (S1)", [self.record("minimap2", "2.24", status="OUTSIDE_ENV",
                                                path="/x/minimap2")])
        with (tracedir / "trace.txt").open("w") as handle:
            handle.write("\t".join(TRACE_HEADER) + "\n")
            for row in self.rows:
                handle.write("\t".join([ "1", row["hash"], "1", row["name"], row["status"], "0",
                                         row["workdir"]]) + "\n")
        for mode, expected in (("strict", 1), ("warn", 0)):
            record = {"invocation": {"params": {"tool_version_check": mode}}, "environment": {}}
            (tracedir / "provenance" / "s_r.json").write_text(json.dumps(record))
            result = subprocess.run(["python3", str(REPO / "scripts" / "finalize_provenance.py"),
                                     "--tracedir", str(tracedir), "--session", "s_r", "--success", "true"],
                                    capture_output=True, text=True, env=clean_env())
            self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
            saved = json.loads((tracedir / "provenance" / "s_r.json").read_text())
            self.assertEqual(saved["environment"]["tool_versions"]["check"], mode)
            report = (tracedir / "RUN_REPORT.txt").read_text()
            self.assertIn("TOOL VERSION MISMATCH", report)
            self.assertIn("[provenance] TOOL VERSION MISMATCH OUTSIDE_ENV", result.stdout)


# ---------------------------------------------------------------------- Nextflow execution

def conda_base():
    conda = shutil.which("conda")
    if not conda:
        return None
    try:
        out = subprocess.run([conda, "info", "--base"], capture_output=True, text=True, timeout=120,
                             env={k: v for k, v in os.environ.items() if k != "LD_LIBRARY_PATH"})
    except (OSError, subprocess.TimeoutExpired):
        return None
    return Path(out.stdout.strip()) if out.returncode == 0 and out.stdout.strip() else None


@unittest.skipUnless(shutil.which("nextflow") and shutil.which("java") and shutil.which("conda"),
                     "needs nextflow, java and conda on PATH")
class NextflowGuard(unittest.TestCase):
    """K06's acceptance test, on Nextflow: a fake conda env pinned to minimap2 2.28, and a fake
    'login' minimap2 2.24 AHEAD of the conda base on the launch PATH (the CRC50 layout)."""

    @classmethod
    def setUpClass(cls):
        cls.base = conda_base()
        if cls.base is None:
            raise unittest.SkipTest("conda info --base failed")
        cls.tmp = Path(tempfile.mkdtemp(prefix="cmp_toolver_nf_"))
        prefix = cls.tmp / "fake_env"
        (prefix / "conda-meta").mkdir(parents=True)
        (prefix / "conda-meta" / "history").write_text("")
        (prefix / "conda-meta" / "minimap2-2.28-h0_0.json").write_text("{}")
        fake_tool(prefix / "bin", "minimap2", "2.28-r1209")
        fake_tool(cls.tmp / "login_tools", "minimap2", "2.24-r1122")
        project = cls.tmp / "project"
        project.mkdir()
        (project / "lib").symlink_to(REPO / "lib")
        (project / "main.nf").write_text(textwrap.dedent(f"""\
            nextflow.enable.dsl=2
            params.scripts = '{REPO / "scripts"}'
            params.tool_version_check = 'strict'
            process probe {{
                conda '{prefix}'
                cpus 1
                memory '200 MB'
                output:
                path 'out.txt'
                script:
                \"\"\"
                ${{ToolVersions.record(task, params, 'minimap2')}}
                minimap2 --version > out.txt
                \"\"\"
            }}
            workflow {{ probe() }}
            """))
        (cls.tmp / "no_guard.config").write_text("process.beforeScript = ''\n")
        cls.project = project
        cls.prefix = prefix

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def run_pipeline(self, name, *extra):
        launch = self.tmp / name
        launch.mkdir()
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("CONDA", "NXF_", "LD_LIBRARY", "PYTHON", "_CE_"))}
        java = Path(shutil.which("java")).parent
        env.update(PATH=":".join([str(self.tmp / "login_tools"), str(self.base / "bin"),
                                  str(self.base / "condabin"), str(java), str(Path(shutil.which("nextflow")).parent),
                                  "/usr/bin", "/bin"]),
                   CONDA_SHLVL="1", CONDA_PREFIX=str(self.base), NXF_ANSI_LOG="false")
        result = subprocess.run(["nextflow", "-q", "run", str(self.project / "main.nf"),
                                 "-c", str(REPO / "nextflow.config"), *extra],
                                cwd=launch, env=env, capture_output=True, text=True, timeout=900)
        records = list(launch.glob("work/*/*/.tool_versions.tsv"))
        rows = []
        if records:
            with records[0].open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        return result, rows

    def test_guard_makes_the_env_win(self):
        result, rows = self.run_pipeline("guarded")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["status"], rows[0]["version"]), ("ok", "2.28"))
        self.assertTrue(rows[0]["path"].startswith(str(self.prefix)), rows[0]["path"])

    def test_without_the_guard_strict_mode_stops_the_run(self):
        result, rows = self.run_pipeline("unguarded", "-c", str(self.tmp / "no_guard.config"))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["status"], rows[0]["version"]), ("OUTSIDE_ENV", "2.24"))
        self.assertIn("78", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
