#!/usr/bin/env bash
# Build the two batch-correction packages that no conda channel provides (r-cqrreg, r-conqur) from
# conda_envs/recipes/ into a new output channel directory. Linux x86_64; needs conda-forge access.
# The committed packages in conda_envs/local_channel/ were built with this script; the batch
# lockfile lists them by checksum, so a rebuild is only needed to change them.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT= CONDA_BIN=${CONDA_EXE:-conda}
while (( $# )); do
  case "$1" in
    --out) OUT=${2:?Missing output directory}; shift 2 ;;
    --conda) CONDA_BIN=${2:?Missing conda executable}; shift 2 ;;
    *) printf 'Usage: bash scripts/build_local_conda_packages.sh --out /new/channel/dir [--conda /path/to/conda]\n' >&2; exit 2 ;;
  esac
done
[[ "$OUT" = /* && ! -e "$OUT" ]] || { printf 'An absolute, not yet existing --out directory is required\n' >&2; exit 2; }
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { printf 'Linux x86_64 only (the batch lockfile platform)\n' >&2; exit 2; }
CONDA_BIN=$(command -v "$CONDA_BIN")
WORK=$(mktemp -d "${TMPDIR:-/tmp}/cmp-conda-build.XXXXXX")
trap 'rm -rf -- "$WORK"' EXIT
# A throwaway builder environment; its package list is kept with the output for provenance.
"$CONDA_BIN" create --yes --quiet --prefix "$WORK/builder" --override-channels -c conda-forge conda-build
"$CONDA_BIN" list --prefix "$WORK/builder" --explicit --md5 > "$WORK/builder.explicit.txt"
BUILD=("$WORK/builder/bin/conda-build" --override-channels -c conda-forge --croot "$WORK/croot"
       --output-folder "$OUT" --no-anaconda-upload)
"${BUILD[@]}" "$ROOT/conda_envs/recipes/r-cqrreg"
# r-conqur needs the r-cqrreg just built: the output folder is an indexed local channel.
"${BUILD[@]}" -c "file://$OUT" "$ROOT/conda_envs/recipes/r-conqur"
cp "$WORK/builder.explicit.txt" "$OUT/builder_environment.explicit.txt"
( cd "$OUT" && find . -name 'r-*.conda' -o -name 'r-*.tar.bz2' | sort | while read -r f; do
    printf '%s  %s  %s\n' "$(md5sum < "$f" | cut -d' ' -f1)" "$(sha256sum < "$f" | cut -d' ' -f1)" "$f"
  done ) | tee "$OUT/packages.md5.sha256"
printf 'LOCAL_CHANNEL_READY=%s\n' "$OUT"
