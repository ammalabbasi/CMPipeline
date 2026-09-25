#!/usr/bin/env bash
# Build an independent Linux x86_64 batch environment from pinned package sources.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PREFIX= CONDA_BIN=${CONDA_EXE:-conda} SOURCE_DIR=
while (( $# )); do
  case "$1" in
    --prefix) PREFIX=${2:?Missing prefix}; shift 2 ;;
    --conda) CONDA_BIN=${2:?Missing conda executable}; shift 2 ;;
    --source-dir) SOURCE_DIR=${2:?Missing source directory}; shift 2 ;;
    *) printf 'Usage: bash scripts/setup_batch_env.sh --prefix /new/env [--conda /path/to/conda] [--source-dir /pinned/tarballs]\n' >&2; exit 2 ;;
  esac
done
[[ "$PREFIX" = /* && "$PREFIX" != / ]] || { printf 'An absolute new --prefix is required\n' >&2; exit 2; }
[[ ! -e "$PREFIX" ]] || { printf 'Refusing to overwrite existing prefix: %s\n' "$PREFIX" >&2; exit 2; }
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { printf 'This lockfile supports Linux x86_64 only\n' >&2; exit 2; }
CONDA_BIN=$(command -v "$CONDA_BIN")
LOCK="$ROOT/conda_envs/lock/batch_correction_env.linux-64.txt"
[[ -s "$LOCK" ]] || { printf 'Missing batch lockfile\n' >&2; exit 2; }
TMP=$(mktemp -d "${TMPDIR:-/tmp}/cmp-batch-sources.XXXXXX")
trap 'rm -rf -- "$TMP"' EXIT
PIN=ff233085cedc36a24382038067a4cece3e94e0a0
CQRTAR=cqrReg_1.2.1.tar.gz
CONQURTAR=ConQuR_par-${PIN}.tar.gz
if [[ -n "$SOURCE_DIR" ]]; then
  cp -- "$SOURCE_DIR/$CQRTAR" "$SOURCE_DIR/$CONQURTAR" "$TMP/"
else
  curl --fail --location --retry 2 "https://cran.r-project.org/src/contrib/Archive/cqrReg/$CQRTAR" -o "$TMP/$CQRTAR"
  curl --fail --location --retry 2 "https://codeload.github.com/ivartb/ConQuR_par/tar.gz/$PIN" -o "$TMP/$CONQURTAR"
fi
cat > "$TMP/SHA256SUMS" <<'SUMS'
7fa9d92c693bfdf155413f4ccbbcc27627f916d290e8bd4944f10c05e82a2647  ConQuR_par-ff233085cedc36a24382038067a4cece3e94e0a0.tar.gz
c235d0ef31e171c7197656bee73429c93c094e80e2ffc856d9b5341caa85089d  cqrReg_1.2.1.tar.gz
SUMS
(cd "$TMP" && sha256sum --check SHA256SUMS)
"$CONDA_BIN" create --yes --prefix "$PREFIX" --file "$LOCK"
# Activate compiler/toolchain variables needed to build the two pinned R packages.
set +u
eval "$("$CONDA_BIN" shell.bash hook)"
conda activate "$PREFIX"
set -u
export R_LIBS="$PREFIX/lib/R/library" R_LIBS_USER="$PREFIX/lib/R/library" R_LIBS_SITE="$PREFIX/lib/R/library"
export CONQUR_REPO=ivartb/ConQuR_par CONQUR_REF="$PIN" CONQUR_OFFLINE_DIR="$TMP"
export MAKEFLAGS=-j1
"$PREFIX/bin/Rscript" --vanilla "$ROOT/scripts/install_conqur.R"
"$PREFIX/bin/Rscript" --vanilla -e 'pkgs <- c("ConQuR", "cqrReg", "vegan", "zCompositions", "doParallel"); x <- data.frame(package=pkgs, version=vapply(pkgs, function(p) as.character(packageVersion(p)), "")); write.table(x, file=file.path(Sys.getenv("CONDA_PREFIX"), "cmp_batch_versions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)'
sha256sum "$LOCK" > "$PREFIX/cmp_batch_lock.sha256"
cp "$TMP/SHA256SUMS" "$PREFIX/cmp_batch_sources.sha256"
printf 'CMP_BATCH_ENV_READY=%s\n' "$PREFIX"
