#!/usr/bin/env python3
"""The batch lockfile's two local packages (r-cqrreg, r-conqur) and conda_envs/local_channel/.

Nextflow builds the batch environment from the lockfile; the two packages no conda channel
provides are listed under @LOCAL_CHANNEL@ and must be the committed files, byte for byte (md5).
The recipes and the module must agree on the pinned ConQuR commit.
"""
import hashlib
import re
import shutil
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "conda_envs" / "lock" / "batch_correction_env.linux-64.txt"
CHANNEL = ROOT / "conda_envs" / "local_channel"
PLACEHOLDER = "@LOCAL_CHANNEL@"


def local_lines():
    return [l for l in LOCK.read_text().splitlines() if l.startswith(PLACEHOLDER)]


class LocalChannel(unittest.TestCase):
    def test_lockfile_ends_with_the_two_local_packages(self):
        lines = LOCK.read_text().rstrip("\n").splitlines()
        local = local_lines()
        self.assertEqual(len(local), 2)
        self.assertEqual(lines[-2:], local, "local packages must come after their dependencies")
        names = [re.search(r"/(r-[a-z]+)-", l).group(1) for l in local]
        self.assertEqual(names, ["r-cqrreg", "r-conqur"])

    def test_listed_files_are_the_committed_ones(self):
        for line in local_lines():
            with self.subTest(line=line):
                rel, _, md5 = line[len(PLACEHOLDER) + 1:].partition("#")
                path = CHANNEL / rel
                self.assertTrue(path.is_file(), f"{rel} is not in conda_envs/local_channel/")
                self.assertRegex(md5, r"^[0-9a-f]{32}$")
                self.assertEqual(hashlib.md5(path.read_bytes()).hexdigest(), md5)

    def test_channel_holds_nothing_unlisted(self):
        listed = {l[len(PLACEHOLDER) + 1:].partition("#")[0] for l in local_lines()}
        present = {p.relative_to(CHANNEL).as_posix() for p in CHANNEL.rglob("*")
                   if p.suffix in (".conda", ".bz2")}
        self.assertEqual(present, listed)

    def test_conqur_pin_is_the_same_everywhere(self):
        main = (ROOT / "main.nf").read_text()
        sha = re.search(r'params\.conqur_sha\s*=\s*"([0-9a-f]{40})"', main).group(1)
        recipe = ROOT / "conda_envs" / "recipes" / "r-conqur"
        self.assertIn(f'{{% set sha = "{sha}" %}}', (recipe / "meta.yaml").read_text())
        self.assertIn(f"CONQUR_SHA={sha}", (recipe / "build.sh").read_text())

    def test_every_resolver_uses_the_same_placeholder(self):
        self.assertIn(f"PLACEHOLDER = '{PLACEHOLDER}'", (ROOT / "lib" / "BatchEnv.groovy").read_text())
        self.assertIn(f"s#{PLACEHOLDER}#file://$ROOT/conda_envs/local_channel#",
                      (ROOT / "scripts" / "setup_batch_env.sh").read_text())
        self.assertIn("params.batch_corr_env = BatchEnv.lockfile(projectDir, workflow.workDir)",
                      (ROOT / "main.nf").read_text())

    def test_same_path_rule_everywhere(self):
        groovy = (ROOT / "lib" / "BatchEnv.groovy").read_text()
        bash = (ROOT / "scripts" / "setup_batch_env.sh").read_text()
        self.assertIn("SAFE_PATH = /[A-Za-z0-9._\\/+-]+/", groovy)
        self.assertIn('[[ "$ROOT" =~ ^[A-Za-z0-9._/+-]+$ ]]', bash)

    def test_vendored_sources_are_the_recipe_sources(self):
        sources = ROOT / "conda_envs" / "recipes" / "sources"
        notice = (CHANNEL / "NOTICE.md").read_text()
        for name in ("r-cqrreg", "r-conqur"):
            meta = (ROOT / "conda_envs" / "recipes" / name / "meta.yaml").read_text()
            sha = re.search(r"sha256: ([0-9a-f]{64})", meta).group(1)
            matches = [p for p in sources.iterdir()
                       if hashlib.sha256(p.read_bytes()).hexdigest() == sha]
            with self.subTest(recipe=name):
                self.assertEqual(len(matches), 1, f"no committed source with the {name} checksum")
                self.assertIn(matches[0].name, notice)
        self.assertEqual(len(list(sources.iterdir())), 2)
        conqur = next(sources.glob("ConQuR_par-*.tar.gz"))
        with tarfile.open(conqur) as tf:
            self.assertIn(f"{conqur.name[:-len('.tar.gz')]}/LICENSE", tf.getnames())
        self.assertIn("GNU GENERAL PUBLIC LICENSE\n                       Version 2",
                      (ROOT / "conda_envs" / "recipes" / "licenses" / "GPL-2").read_text())

    def test_setup_script_refuses_an_unsafe_checkout_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "check#out"
            (root / "scripts").mkdir(parents=True)
            shutil.copy(ROOT / "scripts" / "setup_batch_env.sh", root / "scripts")
            r = subprocess.run(["bash", str(root / "scripts" / "setup_batch_env.sh"),
                                "--prefix", str(Path(tmp) / "env")], capture_output=True, text=True)
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("checkout path may contain only", r.stderr)


@unittest.skipUnless(shutil.which("nextflow"), "nextflow not on PATH")
class ResolvedLockfile(unittest.TestCase):
    """lib/BatchEnv.groovy run by Nextflow: one immutable file per content, unsafe paths refused."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.work = self.tmp / "work"

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def checkout(self, name):
        d = self.tmp / name
        (d / "lib").mkdir(parents=True)
        (d / "conda_envs" / "lock").mkdir(parents=True)
        shutil.copy(ROOT / "lib" / "BatchEnv.groovy", d / "lib")
        shutil.copy(LOCK, d / "conda_envs" / "lock")
        (d / "main.nf").write_text(
            'workflow {\n    println "LOCK=" + BatchEnv.lockfile(projectDir, workflow.workDir)\n}\n')
        return d

    def run_nf(self, d):
        return subprocess.run(["nextflow", "-q", "run", str(d / "main.nf"), "-work-dir", str(self.work)],
                              cwd=self.tmp, capture_output=True, text=True)

    def lockpath(self, r):
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return Path(re.search(r"^LOCK=(.+)$", r.stdout, re.M).group(1))

    def test_checkouts_sharing_a_work_dir_get_their_own_file(self):
        a, b = self.checkout("a"), self.checkout("b")
        pa = self.lockpath(self.run_nf(a))
        mtime = pa.stat().st_mtime_ns
        pb = self.lockpath(self.run_nf(b))
        self.assertEqual(self.lockpath(self.run_nf(a)), pa)
        self.assertEqual(pa.stat().st_mtime_ns, mtime, "an existing resolved lockfile was rewritten")
        self.assertNotEqual(pa, pb)
        self.assertEqual(pa.parent.name, "conda_locks")
        for p, d in ((pa, a), (pb, b)):
            text = p.read_text()
            self.assertNotIn(PLACEHOLDER, text)
            channel = re.search(r"^(file:///\S+)/noarch/r-conqur-", text, re.M).group(1)
            self.assertTrue(channel.endswith(f"/{d.name}/conda_envs/local_channel"), channel)
            digest = hashlib.sha256(text.encode()).hexdigest()[:16]
            self.assertEqual(p.name, f"batch_correction_env.linux-64.{digest}.txt")

    def test_unsafe_checkout_path_stops_the_run(self):
        for name in ("with space", "hash#dir", "amp&dir"):
            with self.subTest(name=name):
                r = self.run_nf(self.checkout(name))
                self.assertNotEqual(r.returncode, 0)
                self.assertIn("checkout path may contain only", r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main()
