nextflow.enable.dsl=2

/*
 * Cohort report: one self-contained HTML page for the run.
 *
 * Reads the QC summary table and the provenance record only -- it does not re-walk the
 * results tree, so the report and the table cannot disagree. No network, no CDN, no
 * JavaScript libraries, so the page still opens years from now from a filesystem.
 */
process cohortReport {

    label 'process_low'
    publishDir("${params.qc_summary_dir}", mode: 'copy')
    conda "${params.consensus_taxa_env}"
    // Re-rendered on every launch (audit C14): it reads THIS launch's provenance record, which a
    // -resume would otherwise skip in favour of the report built from an earlier launch. Cheap.
    cache false

    input:
    path(qc_summary)
    path(validation_report)     // [] unless alignment validation ran
    path(consensus_status)      // [] unless consensus ran (audit A09)

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    path("cmp_cohort_report.html"), emit: report

    script:
    def validation_arg = validation_report ? "--validation-report \"${validation_report}\"" : ""
    def consensus_arg = consensus_status ? "--consensus-status \"${consensus_status}\"" : ""
    // This launch's provenance record, by exact name (<sessionId>_<runName>.json under tracedir),
    // not "newest file in the directory" (audit C14). Outside the work dir, so passed as a path.
    def record = RunRecord.path(params.tracedir, "${workflow.sessionId}_${workflow.runName}")
    """
    set -euo pipefail
    python3 "${params.scripts}/build_cohort_report.py" \\
        --qc-summary "${qc_summary}" \\
        --output cmp_cohort_report.html \\
        --provenance "${record}" ${validation_arg} ${consensus_arg}
    """
}
