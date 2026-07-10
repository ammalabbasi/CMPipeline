nextflow.enable.dsl=2

process DAConcordance {

    tag "${level_model}"
    label 'process_medium'

    conda "${params.da_concordance_env}"
    publishDir "${params.da_concordance_dir}", mode: 'copy', overwrite: true

    input:
    // Each element: [level_model_key, [list of result paths from all enabled tools]]
    tuple val(level_model), path(result_files)

    output:
    path "da_concordance_${level_model}.tsv",         emit: concordance
    path "da_concordance_${level_model}_pairwise.tsv", optional: true, emit: pairwise

    script:
    def level = level_model.split("_")[0]
    def maaslin_flag  = params.run_maaslin  ? "--maaslin  \$(ls *maaslin_full.tsv    2>/dev/null | head -1)" : ""
    def ancombc_flag  = params.run_ancombc  ? "--ancombc  \$(ls *ancombc2_full.tsv   2>/dev/null | head -1)" : ""
    def birdman_flag  = params.run_birdman  ? "--birdman  \$(ls *birdman_results.tsv 2>/dev/null | head -1)" : ""
    def qadabra_flag  = params.run_qadabra  ? "--qadabra  \$(ls concatenated_differentials.tsv 2>/dev/null | head -1)" : ""
    """
    python ${params.da_concordance_script} \\
        ${maaslin_flag} \\
        ${ancombc_flag} \\
        ${birdman_flag} \\
        ${qadabra_flag} \\
        --level ${level} \\
        --output da_concordance_${level_model}.tsv
    """
}
