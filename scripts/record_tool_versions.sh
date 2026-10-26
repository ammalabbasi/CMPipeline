#!/usr/bin/env bash
# Record the version and path of every tool a task is about to run, and refuse to run one
# that is not the pinned build from the task's own conda environment (audit K06/K07/X04;
# ported from pksProfiler F16, commits 2bc7ca9 and 2c0f8c8).
#
# RUN_REPORT's CONDA ENVS section lists what each lockfile SPECIFIES. It cannot see what ran.
# CRC50 attempt 1 (Job 12380606) ran minimap2 2.24 from ~/.bashrc's tools/minimap2 in all 50
# mapReads tasks while the lockfile and the report said 2.28: the launching shell had tool
# directories ahead of the conda base on PATH, and conda activation keeps them there
# (nextflow.config, process.beforeScript). This asks each binary itself, from inside the
# task, after conda activation.
#
#   record_tool_versions.sh MODE ENV PROCESS OUT TOOL[:PACKAGE] ...
#
#   MODE     strict | warn (--tool_version_check)
#   ENV      the task's conda directive (task.conda): a `conda list --explicit` lockfile
#            (conda_envs/lock/*.txt), a .yml, or an existing environment prefix. Pins come
#            from the lockfile's package URLs, the yml's `- package=version` lines, or the
#            prefix's conda-meta. '' when the process has no conda directive.
#   PROCESS  task.process
#   OUT      the record, as TSV: the task work dir's .tool_versions.tsv, which
#            scripts/finalize_provenance.py reads through the trace's workdir column
#   TOOL     a command the task runs; PACKAGE is the conda package that provides it when the
#            names differ (bgzip:htslib, Rscript:r-base, python3:python)
#
# Per tool, one of these statuses:
#   ok                inside $CONDA_PREFIX and the pinned version
#   unpinned          inside $CONDA_PREFIX; the spec pins no exact version for the package
#   version_unread    inside $CONDA_PREFIX; its version output held no version number
#   OUTSIDE_ENV       resolved somewhere other than the task's conda environment (a directory
#                     earlier on PATH, a shell function or alias), or $CONDA_PREFIX is not the
#                     prefix the task asked for
#   VERSION_MISMATCH  reports a version other than its pin
#   NOT_FOUND         not on PATH at all
#   NO_ENV            no conda environment is active in this task
#
# The upper-case statuses stop the task in strict mode with exit 78 (EX_CONFIG), which
# conf/base.config turns into 'terminate', also for per-sample processes under
# --sample_failure_strategy ignore: the environment is run-wide, so every other sample would
# hit the same problem. In warn mode they are printed and recorded, and RUN_REPORT lists them
# under TOOL VERSION MISMATCH.
#
# No `set -e`: a tool whose version flag exits non-zero is recorded, not fatal. Every pipe ends
# in `sed -n` or `grep ... || true`, never `head`, which under pipefail can make the writer die
# of SIGPIPE.
#
# Reading a version is retried (pksProfiler rc1.1, 2c0f8c8). On a 1,304-sample cohort a few
# tasks per tool printed no version line from a single call while every other task of the same
# tool and env read fine: a cold start from a shared-filesystem conda env as the whole cohort
# starts at once outlived a fixed timeout. Each read:
#   - captures stdout and stderr together (fastp and others print to stderr);
#   - runs under `timeout -k 5`, CMP_TOOL_VERSION_TIMEOUT seconds (default 60) on the first
#     attempt and that times the attempt number after;
#   - is retried, up to CMP_TOOL_VERSION_ATTEMPTS (default 3) in all, after a backoff of
#     CMP_TOOL_VERSION_BACKOFF (default 2) times the attempt number seconds, when it exited
#     non-zero or printed nothing and no version was found. A call that exited 0 with output
#     but no version number is not retried: that is deterministic;
#   - records `attempts` and `read_exit` (124 = timed out), and for version_unread puts the exit
#     status and the first output line in `raw`.
# version_unread never fails the task by itself.
set -uo pipefail

if [[ $# -lt 4 ]]; then
    echo "usage: record_tool_versions.sh MODE ENV PROCESS OUT TOOL[:PACKAGE] ..." >&2
    exit 2
fi

MODE=$1
ENV_SPEC=$2
PROCESS=$3
OUT=$4
shift 4

if [[ "$MODE" != "strict" && "$MODE" != "warn" ]]; then
    echo "ERROR: record_tool_versions.sh: MODE must be strict or warn (got '$MODE')" >&2
    exit 2
fi

PREFIX=${CONDA_PREFIX:-}
PREFIX=${PREFIX%/}

# A non-negative integer from the environment, else the default.
setting() {
    local value=${!1:-}
    [[ "$value" =~ ^[0-9]+$ ]] && printf '%s' "$value" || printf '%s' "$2"
}
READ_ATTEMPTS=$(setting CMP_TOOL_VERSION_ATTEMPTS 3)
(( READ_ATTEMPTS >= 1 )) || READ_ATTEMPTS=1
READ_TIMEOUT=$(setting CMP_TOOL_VERSION_TIMEOUT 60)
(( READ_TIMEOUT >= 1 )) || READ_TIMEOUT=1
READ_BACKOFF=$(setting CMP_TOOL_VERSION_BACKOFF 2)
HAVE_TIMEOUT=0
command -v timeout >/dev/null 2>&1 && HAVE_TIMEOUT=1
TIMEOUT=()

# What kind of spec the task's conda directive is.
SPEC_KIND=none
if [[ -n "$ENV_SPEC" && -d "$ENV_SPEC" ]]; then
    SPEC_KIND=prefix
elif [[ -f "$ENV_SPEC" && ( "$ENV_SPEC" == *.yml || "$ENV_SPEC" == *.yaml ) ]]; then
    SPEC_KIND=yml
elif [[ -f "$ENV_SPEC" ]] && grep -q '^@EXPLICIT' "$ENV_SPEC" 2>/dev/null; then
    SPEC_KIND=lockfile
fi

# The exact version the task's spec pins for PACKAGE, or "" when it pins none.
pin_for() {
    local package=$1 spec
    case "$SPEC_KIND" in
        lockfile)
            # .../bioconda/linux-64/minimap2-2.28-h577a1d6_4.tar.bz2#<md5>: name-version-build.
            # Names may contain '-', versions and builds never do.
            awk -v p="$package" '
                /^[[:space:]]*(#|@|$)/ { next }
                {
                    file = $0
                    sub(/#.*$/, "", file)
                    n = split(file, parts, "/")
                    file = parts[n]
                    sub(/\.(conda|tar\.bz2)$/, "", file)
                    m = split(file, f, "-")
                    if (m < 3) next
                    name = f[1]
                    for (i = 2; i <= m - 2; i++) name = name "-" f[i]
                    if (name == p) { print f[m - 1]; exit }
                }' "$ENV_SPEC"
            ;;
        prefix)
            # conda-meta/<name>-<version>-<build>.json: what conda installed in that prefix.
            local meta
            for meta in "$ENV_SPEC"/conda-meta/"$package"-*.json; do
                [[ -f "$meta" ]] || continue
                meta=${meta##*/}
                meta=${meta%.json}
                local rest=${meta#"$package"-}
                # rest = version-build, exactly one '-' left for the right package.
                [[ "$rest" == *-* && "${rest%-*}" != *-* ]] || continue
                printf '%s' "${rest%-*}"
                return 0
            done
            ;;
        yml)
            spec=$(awk -v p="$package" '
                {
                    line = $0
                    sub(/#.*/, "", line)
                    if (!match(line, /^[[:space:]]*-[[:space:]]*/)) next
                    spec = substr(line, RLENGTH + 1)
                    gsub(/[[:space:]]/, "", spec)
                    sub(/^.*::/, "", spec)
                    name = spec
                    sub(/[=<>!~].*$/, "", name)
                    if (name == p) { print substr(spec, length(name) + 1); exit }
                }' "$ENV_SPEC")
            # "=2.28" and "==2.28" are exact; "=2.28=h1234_0" carries a build after the version.
            case "$spec" in
                ==*) spec=${spec#==} ;;
                =*)  spec=${spec#=} ;;
                *)   return 0 ;;
            esac
            spec=${spec%%=*}
            spec=${spec%.\*}
            spec=${spec%\*}
            [[ "$spec" =~ ^[0-9A-Za-z._+-]+$ ]] && printf '%s' "$spec"
            ;;
    esac
    return 0
}

# A few conda builds print a version other than their package version. scripts/
# tool_version_aliases.tsv lists each such build; finalize_provenance.py reads the same file.
ALIASES="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tool_version_aliases.tsv"

# The version PACKAGE's binary is known to print at conda version PIN, or "".
alias_for() {
    [[ -n "$2" && -f "$ALIASES" ]] || return 0
    awk -F '\t' -v p="$1" -v v="$2" '!/^#/ && $1 == p && $2 == v { print $3; exit }' "$ALIASES"
}

# True when VERSION satisfies PIN: the pin, a more specific release of it, or its alias.
satisfies() {
    local version=$1 pin=$2 package=$3 alias
    [[ "$version" == "$pin" || "$version" == "$pin".* ]] && return 0
    alias=$(alias_for "$package" "$pin")
    [[ -n "$alias" && ( "$version" == "$alias" || "$version" == "$alias".* ) ]]
}

# What the tool says about its own version. Most take --version; these do not.
version_output() {
    local tool=$1
    case "$tool" in
        bracken)                        "${TIMEOUT[@]}" bracken -v ;;
        *)                              "${TIMEOUT[@]}" "$tool" --version ;;
    esac
}

clean() { tr '\t\r\n' '   ' | cut -c1-200; }

# Paths first: bowtie2 prints its own absolute path, and an env named python3.10 would
# otherwise read as a version.
strip_paths() { sed -E 's#[^[:space:]]*/[^[:space:]]*##g'; }

# Ask TOOL for its version, retrying a failed call. Sets VERSION, RAW, ATTEMPTS, READ_EXIT.
read_version() {
    local tool=$1 attempt out rc
    VERSION="" RAW="" ATTEMPTS=0 READ_EXIT=""
    for (( attempt = 1; attempt <= READ_ATTEMPTS; attempt++ )); do
        TIMEOUT=()
        (( HAVE_TIMEOUT )) && TIMEOUT=(timeout -k 5 $(( READ_TIMEOUT * attempt )))
        # Both streams; sed reads to the end, so the tool never sees SIGPIPE.
        out=$(version_output "$tool" </dev/null 2>&1 | sed -n '1,200p'; exit "${PIPESTATUS[0]}")
        rc=$?
        ATTEMPTS=$attempt
        READ_EXIT=$rc
        out=$(printf '%s\n' "$out" | strip_paths)
        RAW=$(printf '%s\n' "$out" | grep -E '[0-9]+\.[0-9]+' | sed -n 1p || true)
        VERSION=$(printf '%s\n' "$RAW" | grep -oE '[0-9]+(\.[0-9]+)+' | sed -n 1p || true)
        [[ -n "$VERSION" ]] && return 0
        # Exit 0 with output but no version number: asking again gives the same answer.
        if (( rc == 0 )) && [[ -n "${out//[[:space:]]/}" ]]; then
            break
        fi
        (( attempt < READ_ATTEMPTS )) && sleep $(( READ_BACKOFF * attempt ))
    done
    local first
    first=$(printf '%s\n' "$out" | grep -vE '^[[:space:]]*$' | sed -n 1p || true)
    RAW="exit $READ_EXIT after $ATTEMPTS attempt(s)$( (( READ_EXIT == 124 )) && printf ' (timed out)')${first:+: $first}"
    return 0
}

# A prefix spec names the environment exactly: the task must be running in it.
WRONG_PREFIX=0
if [[ "$SPEC_KIND" == prefix && -n "$PREFIX" ]]; then
    want=$(cd "$ENV_SPEC" 2>/dev/null && pwd -P)
    have=$(cd "$PREFIX" 2>/dev/null && pwd -P)
    [[ -n "$want" && "$want" != "$have" ]] && WRONG_PREFIX=1
fi

printf 'process\ttool\tpackage\tstatus\tversion\tpinned\tpath\tconda_prefix\tenv\traw\tattempts\tread_exit\n' > "$OUT"

problems=()
for entry in "$@"; do
    tool=${entry%%:*}
    package=$tool
    [[ "$entry" == *:* ]] && package=${entry#*:}
    pinned=$(pin_for "$package")

    path=$(command -v "$tool" 2>/dev/null || true)
    # A shell function or alias (exported into this script) runs instead of any binary.
    if [[ -n "$path" && "$path" != /* ]]; then
        path="($(type -t "$tool" 2>/dev/null || echo shell) $tool)"
    fi
    version=""
    raw=""
    attempts=0
    read_exit=""
    if [[ -z "$path" ]]; then
        status=NOT_FOUND
    else
        read_version "$tool"
        version=$VERSION raw=$RAW attempts=$ATTEMPTS read_exit=$READ_EXIT
        if [[ -z "$PREFIX" ]]; then
            status=NO_ENV
        elif [[ "$path" != "$PREFIX"/* ]] || (( WRONG_PREFIX )); then
            status=OUTSIDE_ENV
        elif [[ -n "$pinned" && -n "$version" ]] && ! satisfies "$version" "$pinned" "$package"; then
            status=VERSION_MISMATCH
        elif [[ -z "$version" ]]; then
            status=version_unread
        elif [[ -z "$pinned" ]]; then
            status=unpinned
        else
            status=ok
        fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$PROCESS" "$tool" "$package" "$status" "$version" "$pinned" "$path" "$PREFIX" \
        "$ENV_SPEC" "$(printf '%s' "$raw" | clean)" "$attempts" "$read_exit" >> "$OUT"

    case "$status" in
        NOT_FOUND)        problems+=("$tool: not found on PATH") ;;
        NO_ENV)           problems+=("$tool: ran from $path (${version:-version unread}) with no conda environment active") ;;
        OUTSIDE_ENV)
            if (( WRONG_PREFIX )); then
                problems+=("$tool: the active environment $PREFIX is not the one this task asked for ($ENV_SPEC)")
            else
                problems+=("$tool: ran from $path (${version:-version unread}), outside this task's conda environment $PREFIX${pinned:+; pinned $package=$pinned}")
            fi ;;
        VERSION_MISMATCH) problems+=("$tool: reports $version from $path; $ENV_SPEC pins $package=$pinned") ;;
    esac
done

if (( ${#problems[@]} )); then
    label=ERROR
    [[ "$MODE" == "warn" ]] && label=WARNING
    {
        echo "$label: tool-version check, $PROCESS (--tool_version_check $MODE):"
        printf '  %s\n' "${problems[@]}"
        echo "  Something earlier on PATH shadows the pinned environment, or the environment is not"
        echo "  the one its lockfile describes. Launch Nextflow from a clean environment (docs/hpc.md);"
        echo "  --tool_version_check warn records this and carries on."
    } >&2
    [[ "$MODE" == "strict" ]] && exit 78
fi
exit 0
