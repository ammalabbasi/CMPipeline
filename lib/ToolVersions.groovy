/**
 * The runtime tool-version record (audit K06/K07/X04; ported from pksProfiler F16): one line at
 * the top of every process that runs a tool from its conda environment.
 *
 *     ${ToolVersions.record(task, params, 'minimap2 samtools')}
 *
 * renders a call to scripts/record_tool_versions.sh, which asks each binary on PATH -- inside
 * the task, after conda activation -- for its version and location, writes them to the task work
 * dir's .tool_versions.tsv, and in strict mode stops the task (exit 78, which conf/base.config
 * turns into 'terminate') when one ran from outside the task's conda environment or differs from
 * its pin. RUN_REPORT's CONDA ENVS section lists what the lockfiles specify; CRC50 attempt 1
 * (Job 12380606) ran minimap2 2.24 in every mapReads task while it said 2.28.
 * scripts/finalize_provenance.py reads the records back through the trace's workdir column.
 *
 * The env passed is the task's own conda directive (task.conda), so the pins compared are always
 * the ones the task was built from: a conda_envs/lock/ lockfile, the batch-correction lockfile
 * BatchEnv writes, or a prefix given with --batch_corr_env.
 *
 * The record goes next to .command.sh rather than into the current directory: under `scratch
 * true` the current directory is node-local and discarded, while $0 is still the work dir's
 * .command.sh. It is a dot file so no output glob picks it up. `|| exit $?` because Nextflow
 * runs .command.sh with -u but not -e, and not every script sets -e before this line.
 *
 * task.ext.tool_version_check overrides --tool_version_check for one process. Only the
 * integration test uses it, for Bracken, whose krakenuniq it replaces with a stub on purpose.
 */
class ToolVersions {

    static final List<String> MODES = ['strict', 'warn']

    /** The exit status record_tool_versions.sh uses in strict mode (EX_CONFIG). */
    static final int EXIT_STATUS = 78

    static String record(task, params, String tools) {
        def env = (task.conda ?: '').toString()
        def mode = (task.ext?.tool_version_check ?: params.tool_version_check).toString()
        return "bash \"${params.scripts}/record_tool_versions.sh\" ${mode} " +
               "'${env}' '${task.process}' \"\$(dirname \"\$0\")/.tool_versions.tsv\" ${tools} || exit \$?"
    }

    /** null when valid, else the message main.nf fails with. */
    static String validate(Object mode) {
        def value = mode?.toString()
        return value in MODES ? null :
            "--tool_version_check must be one of ${MODES.join(', ')} (got ${mode})"
    }
}
