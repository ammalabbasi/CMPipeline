#!/usr/bin/env bash
# Optional: pre-build the batch-correction environment into a prefix, for --batch_corr_env.
# Not needed for a normal run: Nextflow builds the same environment from the same lockfile.
# The lockfile's two local packages (r-cqrreg, r-conqur) come from conda_envs/local_channel/.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PREFIX= CONDA_BIN=${CONDA_EXE:-conda}
while (( $# )); do
  case "$1" in
    --prefix) PREFIX=${2:?Missing prefix}; shift 2 ;;
    --conda) CONDA_BIN=${2:?Missing conda executable}; shift 2 ;;
    *) printf 'Usage: bash scripts/setup_batch_env.sh --prefix /new/env [--conda /path/to/conda]\n' >&2; exit 2 ;;
  esac
done
[[ "$PREFIX" = /* && "$PREFIX" != / ]] || { printf 'An absolute new --prefix is required\n' >&2; exit 2; }
[[ ! -e "$PREFIX" ]] || { printf 'Refusing to overwrite existing prefix: %s\n' "$PREFIX" >&2; exit 2; }
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { printf 'This lockfile supports Linux x86_64 only\n' >&2; exit 2; }
CONDA_BIN=$(command -v "$CONDA_BIN")
LOCK="$ROOT/conda_envs/lock/batch_correction_env.linux-64.txt"
grep -q '@LOCAL_CHANNEL@' "$LOCK" || { printf 'Batch lockfile does not reference the local channel\n' >&2; exit 2; }
TMP=$(mktemp -d "${TMPDIR:-/tmp}/cmp-batch-lock.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
# The same substitution lib/BatchEnv.groovy makes for Nextflow.
sed "s#@LOCAL_CHANNEL@#file://$ROOT/conda_envs/local_channel#" "$LOCK" > "$TMP/batch_correction_env.linux-64.txt"
"$CONDA_BIN" create --yes --prefix "$PREFIX" --file "$TMP/batch_correction_env.linux-64.txt"
CMP_PREFIX="$PREFIX" "$PREFIX/bin/Rscript" --vanilla -e 'pkgs <- c("ConQuR", "cqrReg", "vegan", "zCompositions", "doParallel"); x <- data.frame(package=pkgs, version=vapply(pkgs, function(p) as.character(packageVersion(p)), "")); write.table(x, file=file.path(Sys.getenv("CMP_PREFIX"), "cmp_batch_versions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)'
sha256sum "$LOCK" > "$PREFIX/cmp_batch_lock.sha256"
printf 'CMP_BATCH_ENV_READY=%s\n' "$PREFIX"
