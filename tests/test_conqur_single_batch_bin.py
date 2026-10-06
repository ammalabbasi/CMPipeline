#!/usr/bin/env python3
"""K04: Tune_ConQuR must not stop a whole rank when a prevalence bin sits in one batch (2026-10-06).

ConQuR's uncorrected-table benchmark called adonis(tab ~ batchid) without the single-batch guard it
uses for corrected candidates. CRC50 Test 3 lost genus batch correction to it. The pipeline calls a
guarded copy (scripts/batch_contract.R). The behavioural test runs only where the batch env exists.
"""
import glob, os, subprocess, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _batch_rscript():
    if os.environ.get("CMP_TEST_BATCH_RSCRIPT"):
        return os.environ["CMP_TEST_BATCH_RSCRIPT"]
    for env in sorted(glob.glob(str(ROOT / ".conda_cache" / "env-*"))):
        if os.path.isdir(os.path.join(env, "lib", "R", "library", "ConQuR")):
            return os.path.join(env, "bin", "Rscript")
    return None


class GuardWiring(unittest.TestCase):
    def test_script_calls_the_guarded_tune(self):
        text = (ROOT / "scripts" / "2500703_batch_correction_normalization.r").read_text()
        self.assertIn("result <- guarded_tune_conqur()(", text)
        self.assertNotIn("result <- Tune_ConQuR(", text)

    def test_modification_is_recorded(self):
        self.assertIn("audit K04", (ROOT / "conda_envs" / "local_channel" / "NOTICE.md").read_text())


R_CHECK = r'''
suppressWarnings(suppressPackageStartupMessages(library(ConQuR)))
source(file.path("%(root)s", "scripts", "batch_contract.R"))
g <- guarded_tune_conqur()
stopifnot(identical(environment(g), environment(ConQuR::Tune_ConQuR)))
txt <- paste(deparse(g), collapse = "\n")
stopifnot(grepl("else NA", txt, fixed = TRUE))
stopifnot(grepl("all(is.na(as.numeric(R2$batch_R2)))", txt, fixed = TRUE))
# A different ConQuR must stop, not run unguarded.
f <- ConQuR::Tune_ConQuR; body(f) <- quote(NULL)
msg <- tryCatch({guarded_tune_conqur(f); "no error"}, error = conditionMessage)
stopifnot(grepl("Cannot apply the K04 Tune_ConQuR guard", msg, fixed = TRUE))
cat("OK\n")
'''


@unittest.skipUnless(_batch_rscript(), "no batch-correction conda env built yet")
class GuardBehaviour(unittest.TestCase):
    def test_guard_applies_to_the_pinned_conqur_and_refuses_others(self):
        with tempfile.TemporaryDirectory() as d:
            script = Path(d) / "check.R"
            script.write_text(R_CHECK % {"root": ROOT})
            r = subprocess.run([_batch_rscript(), "--vanilla", str(script)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr[-800:])
        self.assertIn("OK", r.stdout)


if __name__ == "__main__":
    unittest.main()
