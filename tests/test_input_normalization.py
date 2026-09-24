#!/usr/bin/env python3
"""Focused, dependency-free checks for normalized workflow interfaces."""

import importlib.util
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "validate_cram_reference", ROOT / "scripts" / "validate_cram_reference.py"
)
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


class InputNormalizationTests(unittest.TestCase):
    def test_modules_use_sample_reads_contract(self):
        for name in ("fastqc.nf", "Bracken.nf", "metaphlan4.nf", "humann3.nf"):
            text = (ROOT / "Modules" / name).read_text()
            self.assertIn("tuple val(sampleID), path(reads)", text, name)
            self.assertNotIn("path(r1_fastq), path(r2_fastq)", text, name)

    def test_main_accepts_single_or_paired_fastq_and_rejects_ambiguity(self):
        text = (ROOT / "main.nf").read_text()
        self.assertIn("fastq_2/r2 is optional", text)
        self.assertIn("has ambiguous inputs", text)

    def test_core_modules_use_task_cpus(self):
        for name in ("extract_reads.nf", "filter_reads.nf", "map_reads.nf",
                     "fastqc.nf", "Bracken.nf", "metaphlan4.nf", "humann3.nf"):
            text = (ROOT / "Modules" / name).read_text()
            self.assertIn("task.cpus", text, name)

    def test_optional_pangenome_falls_through(self):
        text = (ROOT / "Modules" / "map_reads.nf").read_text()
        self.assertIn("params.pangenome_db ?: ''", text)
        self.assertIn('cp "${sampleID}.hg38.t2t.fastq.gz" "${sampleID}.host_depleted.fastq.gz"', text)

    @patch.object(VALIDATOR.subprocess, "run")
    def test_cram_reference_match(self, run):
        run.side_effect = [
            subprocess.CompletedProcess([], 0, "@SQ\tSN:chr1\tLN:10\tM5:abc\n", ""),
            subprocess.CompletedProcess([], 0, "@SQ\tSN:chr1\tLN:10\tM5:abc\n", ""),
        ]
        VALIDATOR.validate("reads.cram", "reference.fa")

    @patch.object(VALIDATOR.subprocess, "run")
    def test_cram_reference_mismatch(self, run):
        run.side_effect = [
            subprocess.CompletedProcess([], 0, "@SQ\tSN:chr1\tLN:10\tM5:abc\n", ""),
            subprocess.CompletedProcess([], 0, "@SQ\tSN:chr1\tLN:11\tM5:def\n", ""),
        ]
        with self.assertRaisesRegex(ValueError, "CRAM/reference mismatch"):
            VALIDATOR.validate("reads.cram", "reference.fa")

    def test_environment_specs_are_minimal(self):
        expected = {
            "samtools_env.yml": {"python=3.10", "samtools=1.22"},
            "fastp_env.yml": {"fastp=0.24.0"},
            "minimap2_env.yml": {"minimap2=2.28", "samtools=1.21"},
        }
        for name, dependencies in expected.items():
            lines = (ROOT / "conda_envs" / name).read_text().splitlines()
            actual = {line.strip()[2:] for line in lines
                      if line.strip().startswith("- ") and "=" in line}
            self.assertEqual(dependencies, actual, name)

    def test_bare_path_flag_is_rejected_before_any_task(self):
        """`--pangenome_db` with no value arrives as Boolean true and used to fail
        inside mapReads, hours into a run."""
        text = (ROOT / "main.nf").read_text()
        validation = text.split("CONDITIONAL WORKFLOW BASED ON ENTRY POINT")[0]
        self.assertIn("expects a path but was given as a bare flag", validation)
        self.assertIn("pathParam('pangenome_db', params.pangenome_db)", validation)
        for name in ("hg38_db", "t2t_phix_db", "kraken_db", "metaphlan_db"):
            self.assertIn(f"requirePathParam('{name}'", validation, name)

    def test_adapter_scan_is_conditional(self):
        module = (ROOT / "Modules" / "filter_reads.nf").read_text()
        self.assertIn("val(trim_adapters)", module)
        self.assertIn("def adapter_args = trim_adapters ?", module)
        self.assertNotIn('fastp -l 45 --adapter_fasta', module)

        text = (ROOT / "main.nf").read_text()
        self.assertIn('params.adapter_trim="auto"', text)
        self.assertIn("adapterTrimMode == 'always'", text)
        self.assertIn("expected auto, always or never", text)

    def test_zero_read_sample_does_not_stop_the_run(self):
        module = (ROOT / "Modules" / "filter_reads.nf").read_text()
        # fastp exits 0 on empty input but leaves a 0-byte file that fails `gzip -t`.
        self.assertIn('if [[ ! -s "\\$FILTERED" ]]; then', module)
        self.assertIn("reads_after_filtering.txt", module)

        text = (ROOT / "main.nf").read_text()
        self.assertIn("readCount.text.trim().toLong() > 0", text)
        self.assertIn("has 0 reads after filtering", text)

    def test_humann3_takes_the_metaphlan_profile_as_an_input(self):
        """It used to read RESULTS/METAPHLAN4 at runtime, racing metaphlan4 itself."""
        module = (ROOT / "Modules" / "humann3.nf").read_text()
        self.assertIn("path(taxonomic_profile)", module)
        self.assertNotIn("params.metaphlan4_dir", module)

        text = (ROOT / "main.nf").read_text()
        # Inner join: HUMAnN3 runs only on samples MetaPhlAn ran on (decided 2026-09-23).
        self.assertIn(".join(METAPHLAN_PROFILES)\n", text)
        self.assertNotIn(".join(METAPHLAN_PROFILES, remainder: true)", text)

    def test_consensus_runs_even_when_no_sample_reached_metaphlan(self):
        """With no MetaPhlAn table, consensus must still run and pass Bracken through,
        not silently skip and leave decontam with no input."""
        text = (ROOT / "main.nf").read_text()
        self.assertIn(".ifEmpty('NO_METAPHLAN')", text)
        self.assertIn("mp == 'NO_METAPHLAN' ? [] : mp", text)
        module = (ROOT / "Modules" / "preprocess_taxa.nf").read_text()
        self.assertIn('path "consensus_status.tsv"', module)
        self.assertIn("def metaphlan_arg = metaphlan_file ?", module)
        script = (ROOT / "scripts" / "compute_consensus_taxa.py").read_text()
        self.assertIn("bracken_passthrough", script)
        self.assertNotIn("sys.exit(f\"ERROR: Consensus requires", script)

    def test_outputs_publish_under_outdir(self):
        text = (ROOT / "main.nf").read_text()
        self.assertIn('params.outdir = "${launchDir}/RESULTS"', text)
        self.assertNotIn('"${projectDir}/RESULTS/', text)

    def test_metaphlan_index_is_a_parameter(self):
        self.assertIn('params.metaphlan_index=', (ROOT / "main.nf").read_text())
        for name in ("metaphlan4.nf", "humann3.nf"):
            module = (ROOT / "Modules" / name).read_text()
            self.assertNotIn("mpa_vJun23", module, name)
            self.assertIn("params.metaphlan_index", module, name)

    def test_config_has_no_top_level_function(self):
        """check_max() is a top-level function, which the Nextflow 25+ config parser
        rejects with "Unexpected input: '('"."""
        def without_comments(text):
            out = []
            for line in text.splitlines():
                stripped = line.strip()
                if stripped.startswith(("//", "*", "/*", "#")):
                    continue
                out.append(line.split("//")[0])
            return "\n".join(out)

        config = (ROOT / "nextflow.config").read_text()
        base = (ROOT / "conf" / "base.config").read_text()
        self.assertNotIn("def check_max", without_comments(config))
        self.assertNotIn("check_max(", without_comments(base))
        self.assertIn("resourceLimits", config)

    def test_conda_cache_is_pinned_outside_the_work_dir(self):
        config = (ROOT / "nextflow.config").read_text()
        self.assertIn("conda.cacheDir", config)

    def test_launch_scripts_are_submittable(self):
        """Slurm stops reading #SBATCH at the first executable line, so directives
        below it are ignored and the job silently gets the 1 h default."""
        for name in ("main.sbatch", "resume_main.sbatch"):
            path = ROOT / "scripts" / "batch_scripts" / name
            lines = path.read_text().splitlines()
            first_executable = next(
                i for i, line in enumerate(lines[1:], start=1)
                if line.strip() and not line.lstrip().startswith("#")
            )
            directives = [i for i, line in enumerate(lines) if line.startswith("#SBATCH")]
            self.assertTrue(directives, f"{name} has no #SBATCH directives")
            self.assertLess(max(directives), first_executable, name)
            for required in ("--account=", "--partition=", "--qos=", "--time="):
                self.assertTrue(any(required in lines[i] for i in directives),
                                f"{name} is missing {required}")
            self.assertNotIn("l1joseph", path.read_text(), name)


if __name__ == "__main__":
    unittest.main()
