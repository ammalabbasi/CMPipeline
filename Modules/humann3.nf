nextflow.enable.dsl=2

process humann3 {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    tag "${sampleID}"
    label 'process_high'
    scratch true
    publishDir "${params.humann3_dir}", mode: 'copy'
    conda "${params.humann3_env}"

    input:
    // `taxonomic_profile` is MetaPhlAn's output for this sample, staged in by the
    // workflow, or an empty list for samples MetaPhlAn was not run on. It used to be
    // read straight out of the published RESULTS/METAPHLAN4 directory, which raced
    // with the metaphlan4 process running on the same channel: whether HUMAnN3 reused
    // a profile or recomputed one depended on task timing.
    tuple val(sampleID), path(reads), path(taxonomic_profile)

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple path("*_genefamilies.tsv"),
          path("*_pathabundance.tsv"),
          path("*_pathcoverage.tsv")

    script:
    // Define the sample name from the input file name
    def SAMPLE_NAME = sampleID
    def read_files = reads instanceof Collection ? reads : [reads]
    def inputs = read_files.collect { "\"${it}\"" }.join(' ')
    def profile_arg = taxonomic_profile ? "--taxonomic-profile \"${taxonomic_profile}\"" : ""

    """
    out1="${SAMPLE_NAME}_genefamilies.tsv"
    out2="${SAMPLE_NAME}_pathabundance.tsv"
    out3="${SAMPLE_NAME}_pathcoverage.tsv"

    # Concatenate and decompress R1 and R2 fastq files in one pass
    gzip -dc ${inputs} > ${SAMPLE_NAME}_combined.fastq


    # Run HUMAnN3
    humann --input ${SAMPLE_NAME}_combined.fastq \\
        --search-mode uniref90 \\
        --nucleotide-database ${params.humann3_nucleotide_db} \\
        --protein-database ${params.humann3_protein_db} \\
        --metaphlan-options "--bowtie2db ${params.metaphlan_db} -x ${params.metaphlan_index}" \\
        --output . \\
        --output-basename ${SAMPLE_NAME} \\
        --threads ${task.cpus} \\
        --input-format fastq \\
        ${profile_arg}

    # Clean up temp directory
    rm -rf ${SAMPLE_NAME}_humann_temp
    """
}


process merge_humann3 {
    label 'process_medium'
    scratch true
    publishDir "${params.humann3_dir}/merged", mode: 'copy'
    conda "${params.humann3_env}"

    input:
    path genefamilies
    path pathabundance
    path pathcoverage

    output:
    path "humann3_genefamilies.tsv"
    path "humann3_pathabundance.tsv"
    path "humann3_pathcoverage.tsv"
    path "humann3_genefamilies_cpm.tsv"
    path "humann3_pathabundance_relab.tsv"

    script:
    """
    # Organize files into separate directories for humann_join_tables
    mkdir -p gf_dir pa_dir pc_dir

    for f in *_genefamilies.tsv; do
        [ -f "\$f" ] && cp "\$f" gf_dir/
    done
    for f in *_pathabundance.tsv; do
        [ -f "\$f" ] && cp "\$f" pa_dir/
    done
    for f in *_pathcoverage.tsv; do
        [ -f "\$f" ] && cp "\$f" pc_dir/
    done

    # Merge per-sample tables into combined tables
    humann_join_tables --input gf_dir --output humann3_genefamilies.tsv --file_name genefamilies
    humann_join_tables --input pa_dir --output humann3_pathabundance.tsv --file_name pathabundance
    humann_join_tables --input pc_dir --output humann3_pathcoverage.tsv --file_name pathcoverage

    # Normalize: gene families to CPM, pathway abundance to relative abundance
    humann_renorm_table --input humann3_genefamilies.tsv --output humann3_genefamilies_cpm.tsv --units cpm
    humann_renorm_table --input humann3_pathabundance.tsv --output humann3_pathabundance_relab.tsv --units relab
    """
}
