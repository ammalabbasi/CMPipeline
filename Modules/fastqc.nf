nextflow.enable.dsl=2

process FASTQC {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    scratch true
    label 'fastqc'
    publishDir("${params.fastqc_dir}", mode: 'copy')
    conda "${params.fastqc_env}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sampleID), path(reads)

    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:
    def read_files = reads instanceof Collection ? reads : [reads]
    def inputs = read_files.collect { "\"${it}\"" }.join(' ')

    """
    fastqc -t ${task.cpus} -o ./ ${inputs}
    """
}
