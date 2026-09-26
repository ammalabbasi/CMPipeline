#!/usr/bin/env python3
"""The batch lockfile's two local packages (r-cqrreg, r-conqur) and conda_envs/local_channel/.

Nextflow builds the batch environment from the lockfile; the two packages no conda channel
provides are listed under @LOCAL_CHANNEL@ and must be the committed files, byte for byte (md5).
The recipes and the module must agree on the pinned ConQuR commit.
"""
import hashlib
import re
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


if __name__ == "__main__":
    unittest.main()
