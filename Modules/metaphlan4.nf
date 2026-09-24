nextflow.enable.dsl=2

process metaphlan4 {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    label 'process_high'
    scratch true
    // Intermediate read/alignment files are published only with --save_intermediates true (audit C11);
    publishDir "${params.metaphlan4_dir}", mode: 'copy', saveAs: { fn -> (params.save_intermediates.toString() == 'true' || !(fn =~ /\.(fastq|fq|fasta|fa|sam|bam|bowtie2\.bz2)(\.gz|\.bz2)?$/)) ? fn : null }
    conda "${params.metaphlan4_env}"  // Set conda environment

    input:
    tuple val(sampleID), path(reads)
 
    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    // sampleID is carried so downstream steps can join this profile back to its
    // sample instead of looking it up by path in the published output directory.
    tuple val(sampleID),
          path("*.metagenome.bowtie2.bz2"),
          path("*.sam.bz2"),
          path("*.profiled_metagenome.txt")


    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = sampleID
    def read_files = reads instanceof Collection ? reads : [reads]
    def inputs = read_files.collect { "\"${it}\"" }.join(' ')

    """
    echo "[DEBUG] PATH is: \$PATH"
    echo "[DEBUG] which metaphlan:"
    which metaphlan || echo "metaphlan not found"
    metaphlan --version || true

    out1="${SAMPLE_NAME}.metagenome.bowtie2.bz2"
    out2="${SAMPLE_NAME}.sam.bz2"
    out3="${SAMPLE_NAME}.profiled_metagenome.txt"

    # Actual execution
    gzip -dc ${inputs} > "${SAMPLE_NAME}.trimmed.fastq"

    # Run MetaPhlAn4
    metaphlan "${SAMPLE_NAME}.trimmed.fastq" \\
        --bowtie2db "${params.metaphlan_db}" \\
        --index ${params.metaphlan_index} \\
        --bowtie2out "${SAMPLE_NAME}.metagenome.bowtie2.bz2" \\
        -s "${SAMPLE_NAME}.sam.bz2" \\
        --nproc ${task.cpus} \\
        --input_type fastq \\
        -o "${SAMPLE_NAME}.profiled_metagenome.txt"
    """
}
