// =============================================================================
// DECONTAMINATION MODULE
// =============================================================================
//
// Statistical identification and removal of contaminants from microbiome data
// using the decontam R package
//
// =============================================================================

nextflow.enable.dsl=2

process Decontamination {

    tag "${prefix}"
    label 'process_medium'

    conda "${params.decontam_env}"
    publishDir "${params.decontam_dir}/${prefix}", mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata)

    output:
    tuple val(prefix),
          path("${prefix}.decontam_pkg_decontaminated.csv"),
          path("decontamination_plots/*"),
          path("${prefix}.decontamination_summary.txt"), emit: decontaminated
    tuple val(prefix),
          path("${prefix}.decontam_pkg_decontaminated.csv"),
          path(metadata), emit: for_batch_correction

    script:
    """
    Rscript ${params.decontam_script} \\
        --otu_table ${otu_table} \\
        --metadata ${metadata} \\
        --prefix ${prefix} \\
        --threshold ${params.decontam_threshold} \\
        --min_prevalence ${params.decontam_min_prevalence} \\
        --min_abundance ${params.decontam_min_abundance} \\
        --min_batches ${params.decontam_min_batches} \\
        --batch_column ${params.batch_column} \\
        --type_column ${params.type_column} \\
        --tumor_values "${params.tumor_values}" \\
        --control_values "${params.control_values}"
    """
}