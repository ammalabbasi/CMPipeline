# Third-party packages in this directory

The two conda packages here are builds of third-party GPL software, redistributed under their
licenses. The complete corresponding source, the build recipes and the license texts are in this
repository.

| Package file | Software | Version | License (DESCRIPTION) | Authors (DESCRIPTION) | Source in this repository |
|---|---|---|---|---|---|
| `linux-64/r-cqrreg-1.2.1-r43_0.conda` | cqrReg (CRAN) | 1.2.1, 2022-06-07 | GPL (>= 2) | Jueyu Gao & Linglong Kong | `conda_envs/recipes/sources/cqrReg_1.2.1.tar.gz` |
| `noarch/r-conqur-2.0-r43_1.conda` | ConQuR, fork `github.com/ivartb/ConQuR_par` | 2.0, commit `ff233085cedc36a24382038067a4cece3e94e0a0` | GPL-3.0 | Wodan Ling and Michael C. Wu | `conda_envs/recipes/sources/ConQuR_par-ff233085cedc36a24382038067a4cece3e94e0a0.tar.gz` |

- The source archives are the upstream files, unchanged; their SHA-256 values are in the recipes
  (`conda_envs/recipes/<name>/meta.yaml`) and checked by `tests/test_local_conda_channel.py`.
  Upstream: `https://cran.r-project.org/src/contrib/cqrReg_1.2.1.tar.gz` (moves to
  `src/contrib/Archive/cqrReg/` when superseded) and
  `https://codeload.github.com/ivartb/ConQuR_par/tar.gz/ff233085cedc36a24382038067a4cece3e94e0a0`.
- License texts: GPL-2 in `conda_envs/recipes/licenses/GPL-2` (copied from R 4.3.3,
  `share/licenses/GPL-2`); GPL-3 is the `LICENSE` file inside the ConQuR source archive.
- Original r43_0 packages were built with `scripts/build_local_conda_packages.sh` from
  `conda_envs/recipes/` on 2026-09-26. The r43_1 ConQuR rebuild is described below.
  The original build did not modify source. At build time (2026-09-26) ConQuR's `DESCRIPTION` got
  `Remote*`/`Github*` fields that record the commit (`conda_envs/recipes/r-conqur/build.sh`);
  cqrReg is built as distributed.
- Both programs come with no warranty; see their licenses.
- **Run-time modification of ConQuR (2026-10-06, audit K04).** `scripts/batch_contract.R`
  (`guarded_tune_conqur`) makes a modified copy of `Tune_ConQuR` in memory each time batch
  correction runs; this guard remains separate from the package compatibility patch below. Two statements are changed. (1) The adonis2
  benchmark of the uncorrected table gets the same single-batch guard that ConQuR already applies
  to its corrected candidates. (2) When no table in a bin can be scored, the uncorrected table is
  kept. The modified function is distributed under GPL-3 as part of this repository.

- **Build-time compatibility patch (2026-11-14, adonis2-terms-v1).** Build r43_1 applies
  `conda_envs/recipes/r-conqur/adonis2-compat.patch` to the unchanged vendored upstream archive.
  Four scorer calls in `Tune_ConQuR` use `vegan::adonis2(method="bray", permutations=999, by="terms")`
  and named R2 extraction, preserving legacy sequential evaluation. `RemoteSha` still identifies
  upstream ff233085; `CMPipelinePatch` separately identifies our GPL-3 modification.
  Other exported ConQuR functions are outside the pipeline's used tuning path and were not migrated.
  The runtime K04 guard matches the patched expression exactly. Batch environment pins vegan2.6-10.
