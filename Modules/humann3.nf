nextflow.enable.dsl=2

process humann3 {
    tag "${SAMPLE_NAME}"
    label 'process_high'
    scratch true
    publishDir "${params.humann3_dir}", mode: 'copy'
    conda "${params.humann3_env}"

    input:
    tuple path(r1_fastq), path(r2_fastq)

    output:
    tuple path("*_genefamilies.tsv"),
          path("*_pathabundance.tsv"),
          path("*_pathcoverage.tsv")

    script:
    // Define the sample name from the input file name
    def SAMPLE_NAME = r1_fastq.baseName.split('\\.')[0]

    """
    out1="${SAMPLE_NAME}_genefamilies.tsv"
    out2="${SAMPLE_NAME}_pathabundance.tsv"
    out3="${SAMPLE_NAME}_pathcoverage.tsv"

    # Skip condition
    if [[ -f "${params.humann3_dir}/\$out1" && -f "${params.humann3_dir}/\$out2" && -f "${params.humann3_dir}/\$out3" ]]; then
        echo "Skipping humann3: Found \$out1, \$out2, \$out3 in publishDir"
        for f in "\$out1" "\$out2" "\$out3"; do
            if [[ ! -f "\$f" ]]; then
                ln -s "${params.humann3_dir}/\$f" . 2>/dev/null || cp "${params.humann3_dir}/\$f" .
            fi
        done
        exit 0
    fi

    # Concatenate and decompress R1 and R2 fastq files in one pass
    zcat ${r1_fastq} ${r2_fastq} > ${SAMPLE_NAME}_combined.fastq

    # Check for existing MetaPhlAn profile (avoids redundant MetaPhlAn run)
    TAXONOMIC_PROFILE_OPT=""
    if [[ -f "${params.metaphlan4_dir}/${SAMPLE_NAME}.profiled_metagenome.txt" ]]; then
        echo "Using existing MetaPhlAn profile from ${params.metaphlan4_dir}"
        TAXONOMIC_PROFILE_OPT="--taxonomic-profile ${params.metaphlan4_dir}/${SAMPLE_NAME}.profiled_metagenome.txt"
    fi

    # Run HUMAnN3
    humann --input ${SAMPLE_NAME}_combined.fastq \\
        --search-mode uniref90 \\
        --nucleotide-database ${params.humann3_nucleotide_db} \\
        --protein-database ${params.humann3_protein_db} \\
        --metaphlan-options "--bowtie2db ${params.metaphlan_db} -x mpa_vJun23_CHOCOPhlAnSGB_202307" \\
        --output . \\
        --output-basename ${SAMPLE_NAME} \\
        --threads ${task.cpus} \\
        --input-format fastq \\
        \$TAXONOMIC_PROFILE_OPT

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
