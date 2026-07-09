nextflow.enable.dsl=2

process MaAsLin {

    tag "${prefix}_${level}_${model}"
    label 'process_high'

    conda "${params.maaslin_env}"
    publishDir "${params.maaslin_dir}/${prefix}_${level}_${model}", mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata), val(level), val(model), val(formula)

    output:
    tuple val("${level}_${model}"), path("*"), emit: results

    script:
    def norm_flag = params.maaslin_normalization ? "--normalization ${params.maaslin_normalization}" : "--normalization TSS"
    def transform_flag = params.maaslin_transform ? "--transform ${params.maaslin_transform}" : "--transform LOG"
    """
    Rscript ${params.maaslin_script} \\
        --otu ${otu_table} \\
        --meta ${metadata} \\
        --prefix ${prefix}_${level} \\
        --formula "${formula}" \\
        --level ${level} \\
        --min_prevalence ${params.da_min_prevalence} \\
        ${norm_flag} \\
        ${transform_flag} \\
        --cores ${task.cpus}
    """
}
