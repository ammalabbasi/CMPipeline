nextflow.enable.dsl=2

process BIRDMAn {

    tag "${prefix}_${level}_${model}"
    label 'process_high'

    conda "${params.birdman_env}"
    publishDir "${params.birdman_dir}/${prefix}_${level}_${model}", mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata), val(level), val(model), val(formula)

    output:
    tuple val("${level}_${model}"), path("*"), emit: results

    script:
    """
    python ${params.birdman_script} \\
        --otu ${otu_table} \\
        --meta ${metadata} \\
        --prefix ${prefix}_${level} \\
        --formula "${formula}" \\
        --level ${level} \\
        --reference ${params.da_reference} \\
        --target ${params.da_target} \\
        --min_prevalence ${params.da_min_prevalence} \\
        --threads ${task.cpus}
    """
}
