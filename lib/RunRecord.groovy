import groovy.json.JsonOutput
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.time.ZonedDateTime
import java.time.format.DateTimeFormatter

/**
 * The equivalence record, launch half: what code, what parameters, what references.
 *
 * Ported from pksProfiler (F16), which had the same gap: nothing tied the released
 * source to the code that actually produced a cohort. Every run now writes one JSON
 * record per session at launch -- before any task runs, so a run that dies still says
 * what it was attempting. The outcome and the output digests are added afterwards by
 * scripts/finalize_provenance.py, called from the workflow.onComplete handler in
 * main.nf.
 *
 * The split is not aesthetic. Classes in lib/ are not visible to the config parser, so
 * the completion half cannot call into here. Keeping it in Python also puts the
 * rendering logic somewhere the test suite can reach.
 *
 * One record per launch, not one per results directory: with `-resume` a results tree
 * is usually the product of several launches, and the case that matters is exactly the
 * one where a later launch ran different code against the same tree. `-resume` reuses
 * the session id, so main.nf keys the record by "<sessionId>_<runName>" (audit C06).
 */
class RunRecord {

    static final int VERSION = 1

    /** Everything knowable before the first task runs. */
    static Map start(Map context) {
        return [
            record_version: VERSION,
            status        : 'started',
            code          : context.code,
            invocation    : context.invocation,
            dependencies  : context.dependencies,
            environment   : context.environment,
            run           : [
                session_id: context.session_id,
                started   : DateTimeFormatter.ISO_OFFSET_DATE_TIME.format(ZonedDateTime.now().withNano(0)),
            ],
            outputs       : [:],
            notes         : context.notes ?: [],
        ]
    }

    /**
     * Where this launch's record lives. `recordKey` is "<sessionId>_<runName>" (main.nf),
     * with no timestamp, so the finaliser can rebuild the same name from `workflow` in
     * onComplete without the two sharing any state. The launch time is inside the record.
     */
    static Path path(Object runsDir, Object recordKey) {
        Path dir = Paths.get(runsDir.toString(), 'provenance')
        Files.createDirectories(dir)
        return dir.resolve("${recordKey}.json")
    }

    static void write(Path target, Map record) {
        Files.createDirectories(target.getParent())
        target.text = JsonOutput.prettyPrint(JsonOutput.toJson(record)) + "\n"
    }

    /**
     * Which stages the parameters asked for. pipeline_info/trace.txt is what actually
     * executed; this is what was requested, which is the part a reader cannot
     * reconstruct afterwards.
     *
     * `flags` carries the already-resolved booleans from the entry workflow rather than
     * re-deriving them here: on Nextflow 26 `--flag false` arrives as the String
     * "false", which is truthy, so there is exactly one place that coercion is done
     * (boolParam in main.nf) and this is not a second one.
     */
    static List<String> stages(Map params, Map flags) {
        def requested = [
            "entry:${params.start_from}".toString(),
            "input:${params.input_data_type}".toString(),
            "adapter_trim:${params.adapter_trim}".toString(),
        ]
        if (params.taxonomic_rank) requested << "rank:${params.taxonomic_rank}".toString()
        if (params.cohort_preset)  requested << "cohort_preset:${params.cohort_preset}".toString()
        if (params.control_role)   requested << "control_role:${params.control_role}".toString()

        // Stages that ran, expressed positively: "skip_humann3 = true" is harder to read
        // in a record than the absence of "humann3".
        // Upstream stages do not run for --start_from decontam / batch_correction.
        def fromStart = params.start_from == 'beginning'
        if (fromStart && !flags.skip_host_depletion) requested << 'host_depletion'
        if (fromStart && !flags.skip_classification) requested << 'classification'
        if (fromStart && !flags.skip_consensus)      requested << 'consensus'
        if (fromStart && !flags.skip_classification && !flags.skip_humann3) requested << 'humann3'
        if (!flags.skip_multiqc)        requested << 'multiqc'
        if (flags.run_decontam && params.start_from != 'batch_correction') requested << 'decontamination'
        if (flags.run_batch_correction) requested << 'batch_correction'
        if (params.pangenome_db)        requested << 'pangenome_depletion'
        if (flags.run_alignment_validation) requested << 'alignment_validation'   // audit C06
        if (!flags.skip_cohort_report)  requested << 'cohort_report'              // audit C06
        return requested
    }

    /**
     * The conda env spec of every stage this launch will run, keyed by parameter name
     * (audit C06/D08). Values are the resolved params.*_env: a lockfile, a yml, or the
     * absolute prefix of a pre-built env (batch correction). The finaliser digests these
     * and reads the prefix env's conda-meta and R DESCRIPTIONs, instead of walking every
     * env in conda.cacheDir -- which lists stale envs and names them by an opaque hash.
     *
     * Mirrors the entry workflow's branches: upstream stages only run from 'beginning',
     * decontamination not from 'batch_correction'. The QC summary and cohort report run
     * in consensus_taxa_env whenever anything feeds them.
     */
    static Map condaEnvs(Map params, Map flags) {
        def names = [] as LinkedHashSet
        boolean fromStart = params.start_from == 'beginning'
        boolean hostDepletion = fromStart && !flags.skip_host_depletion
        boolean classification = fromStart && !flags.skip_classification
        boolean decontam = flags.run_decontam && params.start_from != 'batch_correction'
        if (hostDepletion) {
            names.addAll(['samtools_env', 'fastp_env', 'fastqc_env', 'minimap2_env'])
            if (!flags.skip_multiqc) names << 'multiqc_env'
        }
        if (classification) {
            names.addAll(['krakenuniq_bracken_env', 'metaphlan4_env'])
            if (!flags.skip_consensus) names << 'consensus_taxa_env'
            if (!flags.skip_humann3)   names << 'humann3_env'
            if (flags.run_alignment_validation) names << 'validation_env'
        }
        if (decontam) names << 'decontam_env'
        if (flags.run_batch_correction) names << 'batch_corr_env'
        if ((hostDepletion || classification || decontam)) names << 'consensus_taxa_env'   // qcSummary, cohortReport
        return names.collectEntries { n -> [(n): params[n]?.toString()] }
    }

    /**
     * Which commit this is. workflow.commitId is set only when Nextflow pulled the
     * project itself; most runs are launched from a clone, where it is null and git is
     * the only source. A dirty working tree is recorded rather than ignored -- "commit
     * X" is not true of a run launched with uncommitted edits, and this record is
     * precisely about claims of that kind.
     */
    static Map gitIdentity(Object projectDir, Object commitId, Object revision) {
        def identity = [
            commit        : commitId?.toString(),
            commit_source : commitId ? 'nextflow' : 'unknown',
            revision      : revision?.toString(),
            dirty         : false,
            modified_files: [],
        ]
        try {
            def dir = projectDir.toString()
            def head = capture(['git', '-C', dir, 'rev-parse', 'HEAD'])
            if (head) {
                if (!identity.commit) {
                    identity.commit = head
                    identity.commit_source = 'git'
                }
                else if (identity.commit != head) {
                    identity.commit_source = "nextflow (git HEAD differs: ${head})".toString()
                }
                def status = capture(['git', '-C', dir, 'status', '--porcelain'])
                if (status) {
                    identity.dirty = true
                    identity.modified_files = status.readLines().collect { it.trim() }
                }
            }
        }
        catch (Exception ignored) {
            // not a git work tree, or no git on PATH: keep what we have
        }
        return identity
    }

    private static String capture(List<String> command) {
        def process = new ProcessBuilder(command).start()
        def output = process.inputStream.text.trim()
        return process.waitFor() == 0 ? output : null
    }
}
