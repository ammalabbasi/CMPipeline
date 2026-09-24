"""Repo hygiene checks (audit C21, D10).

Static checks over the pipeline's own source and conda specs. Nothing here reads sample data:
the lockfiles are `conda list --explicit` package-URL lists.

  - every conda_envs/*.yml (top level; conda_envs/unused/ is excluded) is used by the pipeline
  - every lockfile a params.*_env points at exists and has no `defaults` (repo.anaconda.com) URL
  - every withLabel in conf/base.config is used by some process
  - every params.X assigned in an included conf/*.config is declared in main.nf or nextflow.config
  - no *.sbatch has an #SBATCH directive after its first executable line (SLURM ignores those)
"""

import re
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ENVS = REPO / "conda_envs"
LOCK = ENVS / "lock"
BUILD_SBATCH = ENVS / "build_envs.sbatch"

# Only the pipeline's own configs: not work/, RESULTS/ or caches.
NF_SOURCES = [REPO / "main.nf", *sorted((REPO / "Modules").glob("*.nf"))]
CONFIGS = [REPO / "nextflow.config", *sorted((REPO / "conf").glob("*.config"))]

# conf/*.config that configure a DIFFERENT pipeline (run via `-c` from a batch script), so their
# params are not CMPipeline params. mag.config -> nf-core/mag (run_mag.sbatch);
# antismash.config -> run_antismash.nf (run_antismash.sbatch).
FOREIGN_CONFIGS = {"mag.config", "antismash.config"}

PARAM_ENV_RE = re.compile(r"params\.(\w+_env)\s*=\s*[\"']([^\"']+)[\"']")


def _read(p):
    return p.read_text()


def _strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)//.*$", "", text)


def env_params():
    """{param name: path value with ${projectDir} resolved}, from main.nf, Modules and conf."""
    out = {}
    for p in NF_SOURCES + CONFIGS:
        for name, value in PARAM_ENV_RE.findall(_strip_comments(_read(p))):
            out.setdefault(name, value.replace("${projectDir}", str(REPO)))
    return out


def sbatch_prebuilt_ymls(env_values):
    """ymls that build_envs.sbatch creates at a fixed prefix which some params.*_env points at
    (the batch-correction env, audit D04: a lockfile-only build would lack ConQuR)."""
    if not BUILD_SBATCH.exists():
        return set()
    text = _read(BUILD_SBATCH)
    shell_vars = dict(re.findall(r"(?m)^(\w+)=(\S+)\s*$", text))
    shell_vars = {k: v.replace("$REPO", shell_vars.get("REPO", "")) for k, v in shell_vars.items()}
    found = set()
    for var, stem in re.findall(r'conda env create[^\n]*-p "\$(\w+)"[^\n]*conda_envs/(\w+)\.yml', text):
        if shell_vars.get(var) in env_values:
            found.add(stem)
    return found


class CondaEnvHygiene(unittest.TestCase):

    def test_every_top_level_yml_is_used(self):
        """audit D10: a yml nobody references implies an analysis the pipeline does not run."""
        values = set(env_params().values())
        prebuilt = sbatch_prebuilt_ymls(values)
        unused = []
        for yml in sorted(ENVS.glob("*.yml")):
            stem = yml.stem
            used = (str(yml) in values
                    or str(LOCK / f"{stem}.linux-64.txt") in values
                    or stem in prebuilt)
            if not used:
                unused.append(yml.name)
        self.assertEqual(unused, [], "unreferenced conda_envs/*.yml; move to conda_envs/unused/")

    def test_env_params_resolve_to_existing_specs(self):
        missing = [f"params.{k} -> {v}" for k, v in env_params().items()
                   if str(REPO) in v and not Path(v).exists()]
        self.assertEqual(missing, [])

    def test_lockfiles_have_no_defaults_channel(self):
        """audit D07: `defaults` (repo.anaconda.com) is licence-restricted and unreproducible."""
        locks = {Path(v) for v in env_params().values() if v.endswith(".linux-64.txt")}
        self.assertTrue(locks, "no params.*_env points at a lockfile")
        for lf in sorted(locks):
            with self.subTest(lockfile=lf.name):
                self.assertTrue(lf.exists(), f"{lf} missing")
                text = _read(lf)
                self.assertIn("@EXPLICIT", text)
                self.assertNotIn("repo.anaconda.com", text)

    def test_extraction_comment_names_the_locked_samtools(self):
        """audit D11: the category-0 routing claim must name the samtools the lockfile installs."""
        m = re.search(r"/samtools-(\d+(?:\.\d+)+)-", _read(LOCK / "samtools_env.linux-64.txt"))
        self.assertIsNotNone(m, "samtools not found in the samtools_env lockfile")
        self.assertIn(f"samtools {m.group(1)}", _read(REPO / "Modules" / "extract_reads.nf"),
                      "Modules/extract_reads.nf does not cite the locked samtools version; re-check the "
                      "routing on it (cmpipeline/CHECK_ROUND2.sh) and update the comment")

    def test_prebuilt_env_lockfiles_have_no_defaults_channel(self):
        """Lockfiles that build_envs.sbatch writes for pre-built envs (batch correction)."""
        stems = sbatch_prebuilt_ymls(set(env_params().values()))
        self.assertTrue(stems, "build_envs.sbatch no longer builds a pre-built env")
        for stem in sorted(stems):
            lf = LOCK / f"{stem}.linux-64.txt"
            with self.subTest(lockfile=lf.name):
                if not lf.exists():
                    self.skipTest(f"{lf.name} not written yet (conda_envs/build_envs.sbatch "
                                  "has not finished section 3)")
                self.assertNotIn("repo.anaconda.com", _read(lf))


class ConfigHygiene(unittest.TestCase):

    def test_every_base_config_label_is_used(self):
        """audit C21: dead label blocks suggest resources no process receives."""
        base = _strip_comments(_read(REPO / "conf" / "base.config"))
        configured = set(re.findall(r"withLabel\s*:\s*['\"]?(\w+)['\"]?", base))
        self.assertTrue(configured)
        used = set()
        for p in NF_SOURCES:
            used |= set(re.findall(r"(?m)^\s*label\s+['\"](\w+)['\"]", _read(p)))
        self.assertEqual(sorted(configured - used), [], "withLabel blocks no process uses")

    def test_conf_params_are_declared(self):
        """A param set only in a site profile is invisible to users of other profiles, and a
        typo there is silently ignored."""
        declared = set()
        for p in [REPO / "main.nf", REPO / "nextflow.config"]:
            text = _strip_comments(_read(p))
            declared |= set(re.findall(r"(?m)^\s*params\.(\w+)\s*=", text))
            for block in re.findall(r"(?ms)^params\s*\{(.*?)^\}", text):
                declared |= set(re.findall(r"(?m)^\s*(\w+)\s*=", block))
        undeclared = []
        for p in CONFIGS[1:]:
            if p.name in FOREIGN_CONFIGS:
                continue
            text = _strip_comments(_read(p))
            assigned = set(re.findall(r"(?m)^\s*params\.(\w+)\s*=", text))
            for block in re.findall(r"(?ms)^params\s*\{(.*?)^\}", text):
                assigned |= set(re.findall(r"(?m)^\s*(\w+)\s*=", block))
            # A site knob defined AND consumed inside its own profile (e.g. biowulf.config's
            # params.biowulf_lscratch_gb, read by its clusterOptions) is self-contained, not a
            # typo'd or shadowed pipeline param.
            read_here = set(re.findall(r"params\.(\w+)(?!\s*=)", text))
            local = {a for a in assigned if a in read_here
                     and len(re.findall(rf"params\.{a}\b", text)) > 1}
            undeclared += [f"{p.name}: params.{a}" for a in sorted(assigned - declared - local)]
        self.assertEqual(undeclared, [])


class SbatchHygiene(unittest.TestCase):

    def test_no_sbatch_directive_after_first_command(self):
        """SLURM silently ignores #SBATCH lines after the first executable line; the job then
        gets the 1-hour default walltime and times out."""
        scripts = sorted((REPO / "scripts" / "batch_scripts").glob("*.sbatch")) + \
            sorted(ENVS.glob("*.sbatch"))
        self.assertTrue(scripts)
        bad = []
        for s in scripts:
            seen_command = False
            for i, line in enumerate(_read(s).splitlines(), 1):
                stripped = line.strip()
                if stripped.startswith("#SBATCH"):
                    if seen_command:
                        bad.append(f"{s.relative_to(REPO)}:{i}")
                elif stripped and not stripped.startswith("#"):
                    seen_command = True
        self.assertEqual(bad, [], "#SBATCH after the first executable line is ignored")


if __name__ == "__main__":
    unittest.main()
