import java.nio.file.Path

/**
 * The QC floor (pksProfiler COHORT_CHECKS A3; CMPipeline port review #34), the parts main.nf needs.
 *
 * A sample is NOT ASSESSABLE -- never a microbial negative -- when its BAM/CRAM has no FLAG 4
 * record, when its @PG chain shows the unmapped reads were filtered upstream (extractReads decides
 * both and emits `floor`), or when fewer than --min_reads_after_depletion reads survive host
 * depletion. --qc_floor_action exclude (default) stops such a sample; flag lets it carry on. The
 * QC summary (scripts/build_qc_summary.py qc_floor / terminal_state) records it either way.
 *
 * Kept out of main.nf on purpose: Nextflow 24.10's parser stores the entry workflow's source as one
 * Java string, which may not exceed 65,535 characters (tests/test_repo_hygiene.py checks the margin).
 */
class QcFloor {

    static final List<String> ACTIONS = ['exclude', 'flag']

    /** null when both parameters are valid, else the error message. */
    static String problem(Map params) {
        def action = params.qc_floor_action
        if (!(action instanceof CharSequence) || !(action.toString() in ACTIONS))
            return "--qc_floor_action must be exclude or flag (got '${action}')".toString()
        def min = params.min_reads_after_depletion?.toString()
        if (min == null || !min.isLong() || min.toLong() < 0)
            return "--min_reads_after_depletion must be a whole number >= 0 (got '${params.min_reads_after_depletion}'); 0 turns the floor off".toString()
        return null
    }

    /** Whether a sample with extractReads' `floor` continues past extraction. */
    static boolean keepExtracted(Object floor, Object action) {
        return floor?.toString() == 'pass' || action?.toString() == 'flag'
    }

    /** Whether a sample with `n` host-depleted reads continues to classification. */
    static boolean keepDepleted(long n, long min, Object action) {
        return n >= min || action?.toString() == 'flag'
    }

    /** (sample, host-depleted reads) from mapReads' <sample>.depletion_counts.tsv. */
    static List depleted(Path counts) {
        def row = counts.readLines().collect { it.split('\t') }.find { it.size() >= 3 && it[1] == 'host_depleted' }
        if (!row) throw new IllegalStateException("${counts.getFileName()}: no host_depleted row")
        return [row[0], row[2].toLong()]
    }
}
