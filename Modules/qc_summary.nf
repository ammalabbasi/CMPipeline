nextflow.enable.dsl=2

/*
 * Cohort QC summary: one row per sample.
 *
 * Reads only artefacts the pipeline already produces -- fastp JSON, the per-stage
 * depletion counts from mapReads, the KrakenUniq report, the Bracken genus/species
 * reports, the MetaPhlAn profile, and decontamination's sample_validation.tsv. A stage
 * that did not run contributes nothing and its columns read NA.
 *
 * Inputs stage into numbered subdirectories: the same filename can arrive from more
 * than one upstream channel, and a flat staging would collide.
 */
process qcSummary {

    label 'process_low'
    publishDir("${params.qc_summary_dir}", mode: 'copy')
    conda "${params.consensus_taxa_env}"

    input:
    path(qc_inputs, stageAs: "qc_?/*")

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    path("cmp.qc.summary.tsv"), emit: summary

    script:
    """
    set -euo pipefail
    python3 "${params.scripts}/build_qc_summary.py" \\
        --inputs . \\
        --output cmp.qc.summary.tsv \\
        --min-reads "${params.consensus_min_reads}"
    """
}
