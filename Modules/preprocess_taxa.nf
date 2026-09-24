nextflow.enable.dsl=2

process process_metaphlan {

    label 'process_high'
    scratch true
    conda "${params.metaphlan4_env}"
    publishDir "${params.metaphlan4_dir}", mode: 'copy'

    input:
    path profiled_files  // Input directory containing *.profiled_metagenome.txt files

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple path("merged_abundance_table.txt"), 
    path("merged_abundance_table_genus.txt"), 
    path("merged_abundance_table_species.txt"),
    path("merged_abundance_table_SGB.txt")

    script:

    // Prepare space-separated file name lists
    def metaphlan_files = profiled_files.findAll { it.name.endsWith('.profiled_metagenome.txt') }*.getName()
    def metaphlan_str = metaphlan_files.join(' ')


    """
    
    # Run the MetaPhlAn script to merge tables
    find -L . -maxdepth 1 -name '*.profiled_metagenome.txt' -printf '%f\n' | sort > metaphlan_files.txt   # one argument, not one per sample (audit C15)
    python3 "${params.scripts}/merge_metaphlan_tables.py" -l metaphlan_files.txt > merged_abundance_table.txt
    # Sample columns are named by merge_metaphlan_tables.py, which strips only the anchored
    # '.profiled_metagenome.txt' suffix. The unanchored sed passes that used to follow (on every
    # line, taxon names included) are gone (audit C16).


    # Genus level
    (head -n 2 merged_abundance_table.txt && \
     grep -E "g__" merged_abundance_table.txt | grep -v "s__") | \
     sed "s/^.*|//g" > merged_abundance_table_genus.txt

    # Species level
    (head -n 2 merged_abundance_table.txt && \
     grep -E "s__" merged_abundance_table.txt | grep -v "t__") | \
     sed "s/^.*|//g" > merged_abundance_table_species.txt

    # SGB level
    (head -n 2 merged_abundance_table.txt && \
     grep -E "t__" merged_abundance_table.txt) | \
     sed "s/^.*|//g" > merged_abundance_table_SGB.txt

    """
}


process process_bracken {

    label 'process_medium'
    scratch true
    conda "${params.krakenuniq_bracken_env}"
    publishDir "${params.krakenuniq_bracken_dir}", mode: 'copy'

    input:
    path bracken_files

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple path("bracken.genus.mpa.report.txt"),
          path("bracken.species.mpa.report.txt")

    script:
    // Prepare space-separated file name lists for genus and species
    def genus_files = bracken_files.findAll { it.name.endsWith('.G.mpa.krakenreport.txt') }*.getName()
    def species_files = bracken_files.findAll { it.name.endsWith('.S.mpa.krakenreport.txt') }*.getName()

    def genus_str = genus_files.join(' ')
    def species_str = species_files.join(' ')

    """

    # File lists come from find (no per-file argv) and go in as ONE argument, so the merge does not
    # hit ARG_MAX at cohort scale (audit C15). Sorted for a deterministic column order.
    find -L . -maxdepth 1 -name '*.G.mpa.krakenreport.txt' -printf '%f\n' | sort > genus_files.txt
    find -L . -maxdepth 1 -name '*.S.mpa.krakenreport.txt' -printf '%f\n' | sort > species_files.txt

    # Combine genus-level files
    if [ -s genus_files.txt ]; then
        python "${params.scripts}"/combine_mpa.py --input-list genus_files.txt --output bracken.genus.mpa.report.txt
    else
        echo "No genus files found." > bracken.genus.mpa.report.txt
    fi

    # Combine species-level files
    if [ -s species_files.txt ]; then
        python "${params.scripts}"/combine_mpa.py --input-list species_files.txt --output bracken.species.mpa.report.txt
    else
        echo "No species files found." > bracken.species.mpa.report.txt
    fi
    """    
}


process consensus_taxa {

    label 'process_medium'
    scratch true
    publishDir "${params.consensus_taxa_dir}", mode: 'copy'
    conda "${params.consensus_taxa_env}"

    input:
    // metaphlan_file is [] when no sample cleared --consensus_min_reads; the script then
    // passes Bracken through unfiltered and says so in consensus_status.tsv.
    tuple path(metaphlan_file), path(bracken_genus_file), path(bracken_species_file)

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    path "bracken.metaphlan.taxa.prop.pdf"
    path "bracken.metaphlan.common.genus.mpa.report.txt"
    path "bracken.metaphlan.common.species.mpa.report.txt"
    path "consensus_status.tsv"
    path "bracken_library_sizes.tsv"

    script:
    def metaphlan_arg = metaphlan_file ? "--metaphlan ${metaphlan_file}" : ''
    """
    python "${params.scripts}"/compute_consensus_taxa.py ${metaphlan_arg} \\
                                     --bracken_genus ${bracken_genus_file} \\
                                     --bracken_species ${bracken_species_file} \\
                                     --min_reads ${params.consensus_min_reads} \\
                                     --output bracken.metaphlan.taxa.prop.pdf \\
                                     --output_common_genus bracken.metaphlan.common.genus.mpa.report.txt \\
                                     --output_common_species bracken.metaphlan.common.species.mpa.report.txt \\
                                     --output_status consensus_status.tsv \\
                                     --output_library_sizes bracken_library_sizes.tsv
    """
}
