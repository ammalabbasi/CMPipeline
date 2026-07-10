nextflow.enable.dsl=2

// ============================================================================
// WORKFLOW CONTROL PARAMETERS
// ============================================================================

// Pipeline control flags
params.skip_host_depletion = false
params.skip_classification = false
params.skip_consensus = false  // Run consensus taxa (intersect MetaPhlAn + Bracken) before decontam/batch_corr
params.skip_humann3 = true           // Skip HUMAnN3 functional profiling (set false to enable)
params.run_decontam = true           // Enable decontamination step
params.run_batch_correction = false   // disable batch correction step
params.run_diff_abund = false        // Enable differential abundance analysis (Step 6)

// Per-tool DA gates (only used when run_diff_abund = true)
params.run_maaslin  = true
params.run_ancombc  = true
params.run_birdman  = true
params.run_qadabra  = false   // Qadabra requires pre-built Snakemake conda envs
params.run_da_concordance = true

// Entry point for resuming workflow
params.start_from = 'beginning'  // Options: 'beginning', 'decontam', 'batch_correction', 'diff_abund'

// Input files for specific entry points
params.consensus_otu_table = null  // For starting from decontam
params.decontam_otu_table = null   // For starting from batch_correction
params.diff_abund_otu_table = null // For starting from diff_abund (corrected count TSV)
params.metadata_file = null  // Required for decontam/batch_correction: --metadata_file <path>
// ============================================================================
// DECONTAMINATION PARAMETERS
// ============================================================================

params.decontam_threshold = 0.1
params.decontam_min_prevalence = 0.02
params.decontam_min_abundance = 5
params.decontam_min_batches = 2
params.batch_column = "shipment_batch"
params.type_column = "Type"
params.tumor_values = "Tumor"
params.control_values = "Normal"

// ============================================================================
// BATCH CORRECTION PARAMETERS
// ============================================================================

params.batch_corr_covariates = "age_diag,sex,bmi"  // Start conservative
params.tumor_only = false  // have matched tumor/normal, keep both
params.phase = "auto"
params.r2_threshold = 0.25

// ============================================================================
// INPUT/OUTPUT PATHS
// ============================================================================

params.sample = "${projectDir}/samples.csv"

// Output directories
params.unmapped_bam_dir = "${projectDir}/RESULTS/UNMAPPED_BAM"
params.mapped_reads_dir = "${projectDir}/RESULTS/MAPPED_READS"
params.fastqc_dir = "${projectDir}/RESULTS/FASTQC"
params.multiqc_dir = "${projectDir}/RESULTS/MULTIQC"
params.krakenuniq_bracken_dir = "${projectDir}/RESULTS/BRACKEN"
params.metaphlan4_dir = "${projectDir}/RESULTS/METAPHLAN4"
params.humann3_dir = "${projectDir}/RESULTS/HUMANN3"
params.consensus_taxa_dir = "${projectDir}/RESULTS/CONSENSUS_TAXA"
params.decontam_dir = "${projectDir}/RESULTS/04_DECONTAMINATION"
params.batch_corr_dir = "${projectDir}/RESULTS/05_BATCH_CORRECTION"
params.ancombc_dir = "${projectDir}/RESULTS/06_ANCOMBC"
params.maaslin_dir = "${projectDir}/RESULTS/07_MAASLIN"
params.birdman_dir = "${projectDir}/RESULTS/08_BIRDMAN"
params.qadabra_dir = "${projectDir}/RESULTS/09_QADABRA"
params.da_concordance_dir = "${projectDir}/RESULTS/10_DA_CONCORDANCE"
params.antismash_dir = "${projectDir}/RESULTS/ANTISMASH"

// ============================================================================
// DATABASES AND REFERENCE FILES
// ============================================================================

params.hg38_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi"
params.t2t_phix_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi"
params.pangenome_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi"
params.kraken_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023"
params.metaphlan_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/metaphlan"
params.humann3_nucleotide_db='/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/databases/humann3/chocophlan/'
params.humann3_protein_db='/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/databases/humann3/uniref/'
params.adapters="${projectDir}/ref/known_adapters.fna"

// ============================================================================
// CONDA ENVIRONMENT PATHS
// ============================================================================

params.samtools_env = "./conda_envs/samtools_env.yml"
params.fastp_env = "./conda_envs/fastp_env.yml"
params.fastqc_env = "./conda_envs/fastqc_env.yml"
params.minimap2_env = "./conda_envs/minimap2_env.yml"
params.multiqc_env = "./conda_envs/multiqc_env.yml"
params.krakenuniq_bracken_env = "./conda_envs/krakenUniq_bracken_env.yml"
params.consensus_taxa_env = "./conda_envs/consensus_taxa_env.yml"
params.metaphlan4_env = "./conda_envs/metaphlan4_env.yml"
params.humann3_env = "./conda_envs/humann3_env.yml"
params.decontam_env = "./conda_envs/decontam_env.yml"
params.batch_corr_env = "./conda_envs/batch_correction_env.yml"
params.maaslin_env = "./conda_envs/maaslin_env.yml"
params.ancombc_env = "./conda_envs/ancombc_env.yml"
params.birdman_env = "./conda_envs/birdman_env.yml"
params.qadabra_env = "./conda_envs/qadabra_env.yml"
params.da_concordance_env = "./conda_envs/da_concordance_env.yml"
params.antismash_env = "./conda_envs/antismash_env.yml"
params.batch_corr_method = "tune"  // Batch correction method: "tune", "combat", "combat_seq"

params.krakentools_pack ="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/KrakenTools"
params.metaphlan4_pack ="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/MetaPhlAn-4.1.1"

// ============================================================================
// SCRIPT PATHS
// ============================================================================

params.scripts ="${projectDir}/scripts"
params.decontam_script = "${projectDir}/scripts/251006_decontamination_ver2.R"
params.batch_corr_script = "${projectDir}/scripts/2500703_batch_correction_normalization.r"
params.maaslin_script = "${projectDir}/scripts/run_maaslin.R"
params.ancombc_script = "${projectDir}/scripts/run_ancombc.R"
params.birdman_script = "${projectDir}/scripts/run_birdman.py"
params.da_concordance_script = "${projectDir}/scripts/run_da_concordance.py"

// ============================================================================
// DIFFERENTIAL ABUNDANCE PARAMETERS
// ============================================================================

params.da_levels = "genus,species"      // Taxonomic levels to run
params.da_formula = "Type"              // Default formula (single, for clean concordance)
params.da_formula_name = "type"         // Short name for the default formula
params.da_formulas = null               // Optional: semicolon-separated list of "name:formula"
params.da_reference = "Normal"          // Reference level of type_column
params.da_target = "Tumor"             // Target level of type_column
params.da_min_prevalence = 0.10         // Minimum taxon prevalence fraction
params.maaslin_normalization = "TSS"    // MaAsLin normalization (TSS, CLR, NONE, ...)
params.maaslin_transform = "LOG"        // MaAsLin transform (LOG, NONE, ...)
params.da_p_adj_method = "holm"         // ANCOMBC p-value adjustment method
params.da_alpha = 0.05                  // ANCOMBC significance threshold

// ============================================================================
// MODULE IMPORTS
// ============================================================================

// Host depletion modules
include { extractReads } from './Modules/extract_reads.nf'
include { FASTQC as FASTQC1 } from './Modules/fastqc.nf'
include { FASTQC as FASTQC2 } from './Modules/fastqc.nf'
include { FASTQC as FASTQCHG38 } from './Modules/fastqc.nf'
include { FASTQC as FASTQCT2T } from './Modules/fastqc.nf'
include { FASTQC as FASTQCPANGENOME } from './Modules/fastqc.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads } from './Modules/map_reads.nf'
// Taxonomic classification modules
include { Bracken } from './Modules/Bracken.nf'
include { metaphlan4 } from './Modules/metaphlan4.nf'
include { process_metaphlan; process_bracken; consensus_taxa } from './Modules/preprocess_taxa.nf'
include { humann3; merge_humann3 } from './Modules/humann3.nf'

// Preprocessing modules (optional)
include { Decontamination } from './Modules/decontamination.nf'
include { BatchCorrection } from './Modules/batch_correction.nf'

// Differential abundance modules (Step 6)
include { MaAsLin } from './Modules/maaslin.nf'
include { ANCOMBC } from './Modules/ancombc.nf'
include { BIRDMAn } from './Modules/birdman.nf'
include { Qadabra } from './Modules/qadabra.nf'
include { DAConcordance } from './Modules/da_concordance.nf'

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {

    // ============================================================================
    // CONDITIONAL WORKFLOW BASED ON ENTRY POINT
    // ============================================================================

    if (params.start_from == 'beginning') {
        // Full pipeline from BAM files

        // ------------------- STEP1: HOST DEPLETION ---------------------- //

        if (!params.skip_host_depletion) {
            // Read and parse the sample sheet
            sample_sheet = nextflow.Channel.fromPath(params.sample)
                .splitCsv(header: true)
                .map { row ->
                    row.subMap('patient', 'bam') // Extract relevant metadata
                }

            // Extract reads from BAM files
            extractReads(sample_sheet).set { UNMAPPED_READS }

            UNMAPPED_READS.multiMap { sampleID, r1, r2 ->
                path_only: tuple(r1, r2)
                whole: tuple(sampleID, r1, r2)
            }
            .set { UNMAPPED_READS_MULTI }

            // Perform FastQC on the extracted fastq files
            FASTQC1(UNMAPPED_READS_MULTI.path_only)

            // Filter poor quality reads using fastp
            filterReads(UNMAPPED_READS_MULTI.whole).set { FILTERED_UNMAPPED_READS }

            FILTERED_UNMAPPED_READS.multiMap { sampleID, r1, r2 ->
                path_only: tuple(r1, r2)
                whole: tuple(sampleID, r1, r2)
            }
            .set { FILTERED_UNMAPPED_READS_MULTI }

            // Perform FastQC on the filtered fastq files
            FASTQC2(FILTERED_UNMAPPED_READS_MULTI.path_only)

            // gather the list of pangenome .mmi files
            def mmiFiles = []
            def dir = new File("${params.pangenome_db}")
            dir.eachFileRecurse (groovy.io.FileType.FILES) { file ->
                if (file.name.endsWith('.mmi')) {
                    mmiFiles << file
                }
            }

            // Processing for both READSs
            mapReads(FILTERED_UNMAPPED_READS_MULTI.whole, mmiFiles).set { MAPPED_READS }

            MAPPED_READS.multiMap { sampleID, r1Hg38, r1T2T, r1Pan, r2Hg38, r2T2T, r2Pan ->
                Hg38: tuple(r1Hg38, r2Hg38)
                T2T: tuple(r1T2T, r2T2T)
                PAN: tuple(r1Pan, r2Pan)
            }
            .set { MAPPED_READS_MULTI }

            // Perform FastQC on the mapped Reads
            FASTQCHG38(MAPPED_READS_MULTI.Hg38)
            FASTQCT2T(MAPPED_READS_MULTI.T2T)
            FASTQCPANGENOME(MAPPED_READS_MULTI.PAN)
        } else {
            log.info "Skipping host depletion step"
        }

        // ------------------- STEP2: TAXONOMIC CLASSIFICATION ---------------------- //

        if (!params.skip_classification) {
            // BRACKEN
            Bracken(MAPPED_READS_MULTI.PAN).set { BRACKEN_OUT }

            // Collect all kraken reports once all samples are done
            BRACKEN_OUT.map { kraken_report, classified_fasta, unclassified_fasta, bracken_reports, bracken_krakenreports, bracken_mpa_reports ->
                tuple(bracken_mpa_reports)
            }
            .flatten()
            .collect()
            .set {Bracken_mpa_files}

            // Run process_bracken once all reports are available
            process_bracken(Bracken_mpa_files).set { BRACKEN_FILES }

            // METAPHLAN4
            metaphlan4(MAPPED_READS_MULTI.PAN).set { METAPHLAN_OUT }

            // Collect all kraken reports once all samples are done
            METAPHLAN_OUT.map { bowtie2_files, sam_files, profiled_metagenomes ->
                tuple(profiled_metagenomes)
            }
            .flatten()
            .collect()
            .set {metaphlan4_files}

            process_metaphlan(metaphlan4_files).set { METAPHLAN_FILES }

            // HUMANN3 (functional profiling - runs in parallel with classification merging)
            if (!params.skip_humann3) {
                humann3(MAPPED_READS_MULTI.PAN).set { HUMANN3_OUT }

                // Collect per-sample outputs for merging
                HUMANN3_OUT.multiMap { genefamilies, pathabundance, pathcoverage ->
                    genefamilies: genefamilies
                    pathabundance: pathabundance
                    pathcoverage: pathcoverage
                }.set { HUMANN3_MULTI }

                merge_humann3(
                    HUMANN3_MULTI.genefamilies.flatten().collect(),
                    HUMANN3_MULTI.pathabundance.flatten().collect(),
                    HUMANN3_MULTI.pathcoverage.flatten().collect()
                )
            }
        } else {
            log.info "Skipping taxonomic classification step"
        }

        // -------------------  STEP3: CONSENSUS TAXA ---------------------- //

        // Always create bracken_file_tuple when classification is run (needed for both consensus and skip_consensus paths)
        if (!params.skip_classification) {
            BRACKEN_FILES
                .map { bracken_genus_file, bracken_species_file -> tuple(bracken_genus_file, bracken_species_file) }
                .set { bracken_file_tuple }
        }

        if (!params.skip_consensus) {
            // Step 1: Get merged table
            METAPHLAN_FILES
                .map { merged_table, merged_genus, merged_species, merged_SGB -> merged_genus }
                .set { metaphlan4_merged_table }

            // Step 2: Combine for consensus (bracken_file_tuple already defined above)
            metaphlan4_merged_table
                .combine(bracken_file_tuple)
                .set { consensus_input }

            consensus_taxa(consensus_input).set { CONSENSUS_OUTPUT }
        } else {
            log.info "Skipping consensus taxa step"
        }

    } // End of: if (params.start_from == 'beginning')

    // ============================================================================
    // STEP 4: DECONTAMINATION (OPTIONAL)
    // ============================================================================

    if (params.run_decontam && params.start_from != 'batch_correction') {
        // Check if starting from decontamination step
        if (params.start_from == 'decontam') {
            // Load data from parameters for decontamination entry point
            if (!params.consensus_otu_table || !params.metadata_file) {
                error "ERROR: When starting from 'decontam', you must provide:\n" +
                      "  --consensus_otu_table <path>\n" +
                      "  --metadata_file <path>"
            }
            Channel.of(tuple(
                'decontam_run',
                file(params.consensus_otu_table),
                file(params.metadata_file)
            )).set { ch_for_decontam }
        } else if (params.start_from == 'beginning' && !params.skip_consensus) {
            if (!params.metadata_file) {
                error "ERROR: When using --run_decontam, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            CONSENSUS_OUTPUT
                .map {
                    tuple('consensus_run',
                          file("${params.consensus_taxa_dir}/bracken.metaphlan.common.genus.mpa.report.txt"),
                          file(params.metadata_file))
                }
                .set { ch_for_decontam }
        } else if (params.start_from == 'beginning' && params.skip_consensus) {
            if (!params.metadata_file) {
                error "ERROR: When using --run_decontam with --skip_consensus, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            log.info "Skipping consensus taxa - using Bracken genus output directly for decontamination"
            BRACKEN_FILES
                .map { bracken_genus_file, bracken_species_file ->
                    tuple('bracken_run', bracken_genus_file, file(params.metadata_file))
                }
                .set { ch_for_decontam }
        } else {
            log.error "Invalid configuration for decontamination"
        }

        Decontamination(ch_for_decontam)
        Decontamination.out.for_batch_correction.set { ch_for_batch_corr }
    }

    // ============================================================================
    // STEP 5: BATCH CORRECTION (OPTIONAL)
    // ============================================================================

    if (params.run_batch_correction) {
        if (params.start_from == 'batch_correction') {
            // Load data from parameters for batch correction entry point
            if (!params.decontam_otu_table || !params.metadata_file) {
                error "ERROR: When starting from 'batch_correction', you must provide:\n" +
                      "  --decontam_otu_table <path>\n" +
                      "  --metadata_file <path>"
            }
            Channel.of(tuple(
                'batch_corr_run',
                file(params.decontam_otu_table),
                file(params.metadata_file)
            )).set { ch_batch_input }
        } else if (params.run_decontam) {
            // Use output from decontamination
            ch_for_batch_corr.set { ch_batch_input }
        } else if (params.start_from == 'beginning' && !params.skip_consensus) {
            // No decontam, use consensus output directly
            if (!params.metadata_file) {
                error "ERROR: When using --run_batch_correction without decontam, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            CONSENSUS_OUTPUT
                .map {
                    tuple('consensus_run',
                          file("${params.consensus_taxa_dir}/bracken.metaphlan.common.genus.mpa.report.txt"),
                          file(params.metadata_file))
                }
                .set { ch_batch_input }
        } else if (params.start_from == 'beginning' && params.skip_consensus) {
            if (!params.metadata_file) {
                error "ERROR: When using --run_batch_correction with --skip_consensus, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            log.info "Skipping consensus taxa - using Bracken genus output directly for batch correction"
            BRACKEN_FILES
                .map { bracken_genus_file, bracken_species_file ->
                    tuple('bracken_run', bracken_genus_file, file(params.metadata_file))
                }
                .set { ch_batch_input }
        } else {
            log.error "Invalid configuration for batch correction"
        }

        BatchCorrection(ch_batch_input)
    }

    // ============================================================================
    // STEP 6: DIFFERENTIAL ABUNDANCE (OPTIONAL)
    // ============================================================================

    if (params.run_diff_abund) {

        // ---- Build (prefix, otu_table, metadata) input channel ----
        if (params.start_from == 'diff_abund') {
            if (!params.diff_abund_otu_table || !params.metadata_file) {
                error "ERROR: When starting from 'diff_abund', you must provide:\n" +
                      "  --diff_abund_otu_table <path>\n" +
                      "  --metadata_file <path>"
            }
            Channel.of(tuple(
                'diff_abund_run',
                file(params.diff_abund_otu_table),
                file(params.metadata_file)
            )).set { ch_da_input }

        } else if (params.run_batch_correction) {
            // Prefer batch-corrected output; BatchCorrection does not emit metadata so re-attach
            if (!params.metadata_file) {
                error "ERROR: --metadata_file is required for diff_abund when using batch correction"
            }
            BatchCorrection.out.corrected
                .map { prefix, corrected_files, norm_files, pcoa_files, permanova ->
                    def files = corrected_files instanceof List ? corrected_files : [corrected_files]
                    def otu = files.find { it.name == 'ConQuR_tuned.tsv' } ?: files[0]
                    tuple(prefix, otu, file(params.metadata_file))
                }
                .set { ch_da_input }

        } else if (params.run_decontam) {
            ch_for_batch_corr.set { ch_da_input }

        } else if (params.decontam_otu_table) {
            Channel.of(tuple(
                'diff_abund_run',
                file(params.decontam_otu_table),
                file(params.metadata_file)
            )).set { ch_da_input }

        } else {
            error "ERROR: run_diff_abund requires one of: run_batch_correction, run_decontam, " +
                  "start_from=diff_abund, or --decontam_otu_table"
        }

        // ---- Build level × formula channel ----
        def level_formula_combos = []
        params.da_levels.tokenize(',').each { lvl ->
            if (params.da_formulas) {
                // Semicolon-separated "name:formula" pairs e.g. "type:Type;type_country:Type+country"
                params.da_formulas.tokenize(';').each { entry ->
                    def parts = entry.trim().tokenize(':')
                    level_formula_combos << [lvl.trim(), parts[0].trim(), parts[1..-1].join(':').trim()]
                }
            } else {
                level_formula_combos << [lvl.trim(), params.da_formula_name, params.da_formula]
            }
        }
        Channel.fromList(level_formula_combos).set { ch_levels_formulas }

        // Combine: (prefix, otu, meta) × (level, model_name, formula)
        ch_da_input
            .combine(ch_levels_formulas)
            .set { ch_da_combined }
        // ch_da_combined: (prefix, otu, meta, level, model_name, formula)

        // ---- Dispatch to enabled tools ----
        def ch_concordance_input = Channel.empty()

        if (params.run_maaslin) {
            MaAsLin(ch_da_combined)
            ch_concordance_input = ch_concordance_input.mix(MaAsLin.out.results)
        }
        if (params.run_ancombc) {
            ANCOMBC(ch_da_combined)
            ch_concordance_input = ch_concordance_input.mix(ANCOMBC.out.results)
        }
        if (params.run_birdman) {
            BIRDMAn(ch_da_combined)
            ch_concordance_input = ch_concordance_input.mix(BIRDMAn.out.results)
        }
        if (params.run_qadabra) {
            Qadabra(ch_da_combined)
            ch_concordance_input = ch_concordance_input.mix(Qadabra.out.results)
        }

        // ---- Concordance: group all tool outputs by level_model key ----
        if (params.run_da_concordance) {
            ch_concordance_input
                .groupTuple()
                .set { ch_concordance_grouped }
            DAConcordance(ch_concordance_grouped)
        }
    }

}

// ============================================================================
// WORKFLOW COMPLETION HANDLERS
// ============================================================================

workflow.onComplete {
    log.info """
    ================================================================================
    Pipeline execution summary
    ================================================================================
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    ================================================================================
    """.stripIndent()
}

workflow.onError {
    log.info """
    ================================================================================
    Pipeline execution error
    ================================================================================
    ${workflow.errorMessage}
    ================================================================================
    """.stripIndent()
}
