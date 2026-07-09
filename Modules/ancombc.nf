nextflow.enable.dsl=2

process ANCOMBC {

    tag "${prefix}_${level}_${model}"
    label 'process_high'

    conda "${params.ancombc_env}"
    publishDir "${params.ancombc_dir}/${prefix}_${level}_${model}", mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata), val(level), val(model), val(formula)

    output:
    tuple val("${level}_${model}"), path("*"), emit: results

    script:
    """
    Rscript ${params.ancombc_script} \\
        --otu ${otu_table} \\
        --meta ${metadata} \\
        --prefix ${prefix}_${level} \\
        --formula "${formula}" \\
        --level ${level} \\
        --p_adj_method ${params.da_p_adj_method ?: 'holm'} \\
        --alpha ${params.da_alpha ?: 0.05} \\
        --prv_cut ${params.da_min_prevalence}
    """
}
