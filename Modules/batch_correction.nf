// =============================================================================
// BATCH CORRECTION MODULE
// =============================================================================
//
// Batch effect correction and normalization of microbiome data using ConQuR
//
// =============================================================================

nextflow.enable.dsl=2

process BatchCorrection {

    tag "${prefix}"
    label 'process_high'

    conda "${params.batch_corr_env}"
    publishDir "${params.batch_corr_dir}/${prefix}", mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata)

    output:
    tuple val(prefix),
          path("corrected/*"),
          path("normalized/*.tsv"),
          path("pcoa_plots/*"),
          path("permanova_summary.tsv"), emit: corrected

    script:
    def tumor_only_flag = params.tumor_only ? "--tumor_only" : ""
    def method_flag = params.batch_corr_method ?: "tune"
    """
    Rscript ${params.batch_corr_script} \\
        --otu ${otu_table} \\
        --meta ${metadata} \\
        --prefix ${prefix} \\
        --batch_column ${params.batch_column} \\
        --covariates ${params.batch_corr_covariates} \\
        --type_column ${params.type_column} \\
        --tumor_value ${params.tumor_values} \\
        ${tumor_only_flag} \\
        --phase ${params.phase} \\
        --r2_threshold ${params.r2_threshold} \\
        --method ${method_flag}
    """
}
