nextflow.enable.dsl=2

/*
 * Run antiSMASH to identify biosynthetic gene clusters (BGCs).
 *
 * Supports two input types:
 *   - Metagenomic assembly contigs (uses prodigal-m for gene finding)
 *   - MAG bins (uses prodigal for single-genome gene finding)
 *
 * Input:  tuple(sample_id, fasta_file)
 * Output: antiSMASH results (HTML, JSON, GBK) per sample/MAG
 */

process antiSMASH {
    label 'process_high'
    scratch true
    publishDir "${params.antismash_dir}/${sample_id}", mode: 'copy'
    conda "${params.antismash_env}"

    // Ignore failures from empty/tiny assemblies with no genes
    errorStrategy { task.exitStatus in [141, 143, 137, 104, 134, 139] ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}/*.json"),     emit: json,     optional: true
    tuple val(sample_id), path("${sample_id}/*.gbk"),      emit: gbk,      optional: true
    tuple val(sample_id), path("${sample_id}/index.html"),  emit: html,     optional: true
    tuple val(sample_id), path("${sample_id}/*"),           emit: all

    script:
    // Use prodigal-m (metagenomic mode) for contigs, prodigal for MAGs
    def genefinder = params.mode == "contigs" ? "prodigal-m" : "prodigal"
    """
    # Skip check: if results already exist in publishDir
    if [[ -f "${params.antismash_dir}/${sample_id}/index.html" ]]; then
        echo "Skipping ${sample_id}: results already exist"
        mkdir -p ${sample_id}
        ln -s ${params.antismash_dir}/${sample_id}/* ${sample_id}/ 2>/dev/null || true
        exit 0
    fi

    # Decompress if gzipped
    INPUT_FASTA="${fasta}"
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c "${fasta}" > "${sample_id}.fa"
        INPUT_FASTA="${sample_id}.fa"
    fi

    # Run antiSMASH
    antismash \\
        --taxon bacteria \\
        --output-dir ${sample_id} \\
        --output-basename ${sample_id} \\
        --cpus ${task.cpus} \\
        --genefinding-tool ${genefinder} \\
        --cb-general \\
        --cb-knownclusters \\
        --cb-subclusters \\
        --asf \\
        --pfam2go \\
        --cc-mibig \\
        "\${INPUT_FASTA}"
    """
}
