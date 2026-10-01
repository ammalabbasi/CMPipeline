nextflow.enable.dsl=2

/*
 * Per-stage survival of the reads the low-score rescue added (--extract_low_score_mapped; Ludmil
 * 1 Oct, L8). Runs only for BAM/CRAM samples whose rescue was on.
 *
 * extractReads writes the rescued reads' names exactly as `samtools fastq -N` wrote them. Names
 * survive the rest of the chain unchanged: fastp keeps them, and mapReads's minimap2 -ax sr +
 * `samtools fastq` restores the /1 /2 on any mates it pairs (READ1/READ2 flags). So the rescued reads
 * are counted by name in the host-depleted FASTQ and in KrakenUniq's classified FASTA. Counts only:
 * one TSV row, no name or sequence is written out.
 */
process lowAsSurvival {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    label 'process_low'
    tag "$sampleID"
    publishDir("${params.qc_summary_dir}/low_as_survival", mode: 'copy')
    conda "${params.consensus_taxa_env}"

    input:
    // classified: KrakenUniq's classified FASTA, or an empty list when classification is skipped.
    tuple val(sampleID), path(low_as_ids), path(host_depleted), path(classified)

    // Digests of the helper scripts this task uses (Provenance.code), as in the other modules (audit C03).
    val deps
    output:
    path("${sampleID}.low_as_survival.tsv"), emit: counts

    script:
    def classified_arg = classified ? "--classified \"${classified}\"" : ""
    """
    set -euo pipefail
    python3 "${params.scripts}/low_score_rescue.py" survival \\
        --ids "${low_as_ids}" --depleted "${host_depleted}" ${classified_arg} \\
        --sample "${sampleID}" --output "${sampleID}.low_as_survival.tsv"
    """
}
