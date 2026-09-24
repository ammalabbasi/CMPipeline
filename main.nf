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
params.skip_multiqc = false          // Cohort QC report from fastp JSON + FastQC
params.skip_cohort_report = false    // Self-contained HTML cohort report
// finish | ignore: what a per-sample failure does to the run (audit C10). finish (default) stops it;
// ignore lets the other samples complete, with the failures listed in the run report.
params.sample_failure_strategy = 'finish'
params.save_intermediates = false   // publish intermediate read files too (default: host-depleted reads only; audit C11)

// Entry point for resuming workflow
params.start_from = 'beginning'  // Options: 'beginning', 'decontam', 'batch_correction'

// Input files for specific entry points
params.consensus_otu_table = null  // For starting from decontam
params.decontam_otu_table = null   // For starting from batch_correction
params.metadata_file = null  // Required for decontam/batch_correction: --metadata_file <path>
// ============================================================================
// DECONTAMINATION PARAMETERS
// ============================================================================

params.decontam_threshold = 0.1
// 0.05 = keep taxa detected in >=5% of samples, applied after the >=5-read abundance
// filter. This matches the manuscript's stated CRC pre-filter (">=5 reads, detected
// in >=5% of samples") and June-young's ver4 default. The pipeline had carried 0.02 since
// the committed base, which was the odd one out: --decontam_min_abundance (5) and
// --decontam_threshold (0.1) already matched the published methods.
// Adopted 2026-09-22 (Ammal). See MERGE_CHECKLIST 3.1.
params.decontam_min_prevalence = 0.05
params.decontam_min_abundance = 5
params.decontam_min_batches = 2
// keep | drop -- a sample whose counts are ALL zero in the decontam input table. Such a sample is
// evidence of absence for the prevalence test, so it is kept as long as its library was not
// empty upstream (per-sample Bracken totals from consensus_taxa, or --decontam_library_sizes).
// Only truly empty libraries are excluded, and named. `drop` restores the ver4 behaviour of
// silently excluding every all-zero sample, which on 2026-09-23 removed all three batch-1
// normals from the balanced CRC run (6 consensus genera, none present in those normals) and so
// removed that batch's controls. Decided 2026-09-23 (Ammal).
params.decontam_zero_total_policy = "keep"
params.decontam_library_sizes = null   // for --start_from decontam: bracken_library_sizes.tsv
// skip | error -- a taxon present in a batch but in < 2 of its samples cannot be tested by the
// prevalence method (decontam returns NA). skip: record it as unevaluable, not flagged in that
// batch; the >= --decontam_min_batches rule then decides from the batches that could test it.
// error: stop (ver4 behaviour). Decided 2026-09-23 (Ammal), after the balanced CRC run hit it.
params.decontam_unevaluable_policy = "skip"
// Minimum Bracken library for an all-zero sample to count as evidence of absence; below it the
// sample is excluded as insufficient_library (audit A06/A18; decided 2026-09-23, Ammal: 1000).
params.decontam_min_library = 1000
// A zero-total sample with unknown library size stops decontam unless this is set (audit A07).
params.decontam_allow_unknown_library = false
// Per-batch eligibility floors (ver4).
params.min_total_per_batch = 5
params.min_positive_per_batch = 1
params.min_control_per_batch = 1
params.batch_column = "shipment_batch"
// Batch correction only, and only when --tumor_only is set: which sample types count as
// tumours for that subsetting. NOT the same question as decontamination's
// --positive_values, which means "not a control" and in TCGA includes Solid Tissue Normal.
// Renamed 2026-09-22 so the two cannot be confused. See MERGE_CHECKLIST 3.2.
params.batch_corr_tumor_values = "Tumor"
// A failed batch correction stops the run (audit A01). true = continue with the raw counts,
// published as corrected/UNCORRECTED_passthrough.tsv, never as ConQuR_tuned.tsv.
params.batch_corr_allow_uncorrected = false
// The ConQuR commit the batch-correction script was written against (ivartb/ConQuR_par fork);
// the module refuses any other install (audit A14).
params.conqur_sha = "ff233085cedc36a24382038067a4cece3e94e0a0"

// --- Decontamination sample-role contract (ver4) -----------------------------------
// The taxonomy column in the OTU table; every other column is treated as a sample.
params.taxon_column = "clade_name"
// genus | species -- selects which consensus/Bracken table feeds decontamination.
// Required whenever decontamination or batch correction runs; no safe default.
params.taxonomic_rank = null
params.sample_id_column = "sampleid"
// error | exclude -- what to do with a sample whose type was not declared.
params.unselected_type_policy = "error"
// null | tcga -- a preset for the three role parameters below.
params.cohort_preset = null
// These three are required when decontamination runs and have no cohort-neutral default:
// naming the positives and the negative controls is a scientific choice, not a default.
// A value given on the command line wins (Nextflow keeps the first assignment), so these
// only apply when the flag was not passed.
// The metadata column holding each sample's role. Neutral default, overridable per
// cohort; --cohort_preset tcga supplies the TCGA column name. A CLI value wins over
// both (Nextflow keeps the first assignment). Decided 2026-09-22 (Ammal); see
// MERGE_CHECKLIST 3.3. A column that exists but means something else is caught by
// --unselected_type_policy error, since its values would not match the declared
// --positive_values / --control_values.
params.type_column = params.cohort_preset == 'tcga' ? 'Sample Type_x' : 'sample_type'
params.control_role    = params.cohort_preset == 'tcga' ? 'surrogate_biological' : null
params.positive_values = params.cohort_preset == 'tcga' ? 'Primary Tumor,Solid Tissue Normal' : null
params.control_values  = params.cohort_preset == 'tcga' ? 'Blood Derived Normal' : null

// ============================================================================
// BATCH CORRECTION PARAMETERS
// ============================================================================

params.batch_corr_covariates = "age_diag,sex,bmi"  // Metadata columns ConQuR adjusts for. At least one is
                                                   // required ("none" stops at launch, decision 16), and each
                                                   // named column must exist in --metadata_file
params.tumor_only = false  // have matched tumor/normal, keep both
params.phase = "auto"
params.r2_threshold = 0.25
params.bracken_read_length = 50  // Must match the read length used to build the Bracken database
// Bracken -t: the minimum reads a taxon needs to be re-estimated -- a THRESHOLD, not a thread
// count. The committed value was 2; an edit had set it to ${task.cpus}, so abundance output
// changed with the CPU count and the retry attempt (audit A04, 2026-09-23).
params.bracken_threshold = 2
params.consensus_min_reads = 100000  // Per-sample minimum KrakenUniq root (classified) reads for MetaPhlAn

// ============================================================================
// INPUT/OUTPUT PATHS
// ============================================================================

params.sample = "${projectDir}/samples.csv"
params.input_data_type = "auto"  // auto | bam | cram | fastq
params.cram_reference = null      // optional matching FASTA for CRAM decoding

// Output directories. Everything publishes under --outdir, which defaults to RESULTS/ in
// the directory the run was launched from, next to that run's work/ and pipeline_info/.
// Until 2026-09-23 these were ${projectDir}/RESULTS/..., so every run from any launch dir
// wrote into (and could overwrite) the same folder inside the repo, and a stage-2 run
// looking under its own launch dir never found stage 1's consensus table.
// A value given on the command line wins (Nextflow keeps the first assignment).
params.outdir = "${launchDir}/RESULTS"
params.unmapped_bam_dir = "${params.outdir}/UNMAPPED_BAM"
params.mapped_reads_dir = "${params.outdir}/MAPPED_READS"
params.fastqc_dir = "${params.outdir}/FASTQC"
params.multiqc_dir = "${params.outdir}/MULTIQC"
params.qc_summary_dir = "${params.outdir}/QC_SUMMARY"
// Mirrors conda.cacheDir in nextflow.config. Declared here so the completion handler
// can read it as a param: on Nextflow 24.10 a `def` inside the entry workflow that
// references System or projectDir fails to compile ("Variable already defined in the
// process scope").
params.conda_cache_dir = System.getenv('NXF_CONDA_CACHEDIR') ?: "${projectDir}/.conda_cache"
params.krakenuniq_bracken_dir = "${params.outdir}/BRACKEN"
params.metaphlan4_dir = "${params.outdir}/METAPHLAN4"
params.humann3_dir = "${params.outdir}/HUMANN3"
params.consensus_taxa_dir = "${params.outdir}/CONSENSUS_TAXA"
params.decontam_dir = "${params.outdir}/04_DECONTAMINATION"
params.batch_corr_dir = "${params.outdir}/05_BATCH_CORRECTION"
params.antismash_dir = "${params.outdir}/ANTISMASH"
params.validation_dir = "${params.outdir}/06_VALIDATION"

// ============================================================================
// ALIGNMENT VALIDATION (opt-in). Design: DESIGN_alignment_validation.md (v2, 2026-09-23).
// Thresholds are the simulation benchmark's (alignment_benchmark/RESULTS.md); identity is the
// least certain and should be revisited after the first real run.
// ============================================================================
params.run_alignment_validation = false
params.validation_min_support = 25        // reads on the best genome; below this: insufficient
params.validation_min_confirmed = 10      // Bracken >= min_support but < this many supporting reads: unsupported, zeroed (decision 29)
params.validation_min_bins_ratio = 0.5    // 10 kb windows hit / expected (per read pair)
params.validation_min_identity = 0.98     // median 1 - NM/aligned columns
params.validation_min_unique_frac = 0.02  // unique / (unique + shared) reads
params.validation_max_genomes = 50        // RefSeq genomes per species, before de-duplication
params.validation_dedup_ani = 99.5        // skani ANI at which two genomes of a species are one
params.validation_max_missing_frac = 0.10 // fail if more species than this get no genome
params.validation_max_ref_gb = 40         // refuse a planned reference above this many Gb, before downloading (audit B05)
params.validation_batch_size = 20         // samples per alignment task (index loads once each)
params.validation_groups = "${projectDir}/ref/validation_groups.tsv"   // E. coli + Shigella
params.validation_export_reads = false    // write reads supporting validated taxa (protected data)
params.validation_taxdb = null            // default: <kraken_db>/taxDB
params.validation_cache_dir = System.getenv('CMP_VALIDATION_CACHE') ?: "${projectDir}/.validation_cache"
// A pre-built reference directory (a ref-<key>/ folder from the cache) skips the species list and
// the build: for runs without internet, or to pin one reference across cohorts.
params.validation_ref = null
// RefSeq snapshot month (YYYY-MM) of the validation reference; part of its cache key. Default: the
// month the run starts. A past month is reused from the cache if there, else the build stops.
params.validation_refseq_snapshot = null

// ============================================================================
// DATABASES AND REFERENCE FILES
// ============================================================================

params.hg38_db=null       // Set in a site profile or with --hg38_db
params.t2t_phix_db=null    // Set in a site profile or with --t2t_phix_db
params.pangenome_db=null  // optional .mmi file or directory containing .mmi files
params.kraken_db=null     // Set in a site profile or with --kraken_db
params.metaphlan_db=null  // Set in a site profile or with --metaphlan_db
params.metaphlan_index="mpa_vJun23_CHOCOPhlAnSGB_202307"  // Bowtie2 index name inside --metaphlan_db
params.humann3_nucleotide_db=null  // Set in a site profile or with --humann3_nucleotide_db
params.humann3_protein_db=null     // Set in a site profile or with --humann3_protein_db
params.adapters="${projectDir}/ref/known_adapters.fna"
params.adapter_trim="auto"  // auto | always | never -- see the adapter_trim note in the workflow below

// ============================================================================
// CONDA ENVIRONMENT PATHS
// ============================================================================

// Exact, locked environments (`conda list --explicit`), one per module (audit D01/D02, decision 23).
// The .yml beside each lockfile is its human-readable source; conda_envs/build_envs.sbatch
// regenerates the lockfiles. No `defaults` channel anywhere (audit D07).
params.samtools_env = "${projectDir}/conda_envs/lock/samtools_env.linux-64.txt"
params.fastp_env = "${projectDir}/conda_envs/lock/fastp_env.linux-64.txt"
params.fastqc_env = "${projectDir}/conda_envs/lock/fastqc_env.linux-64.txt"
params.minimap2_env = "${projectDir}/conda_envs/lock/minimap2_env.linux-64.txt"
params.multiqc_env = "${projectDir}/conda_envs/lock/multiqc_env.linux-64.txt"
params.krakenuniq_bracken_env = "${projectDir}/conda_envs/lock/krakenUniq_bracken_env.linux-64.txt"
params.consensus_taxa_env = "${projectDir}/conda_envs/lock/consensus_taxa_env.linux-64.txt"
params.metaphlan4_env = "${projectDir}/conda_envs/lock/metaphlan4_env.linux-64.txt"
params.humann3_env = "${projectDir}/conda_envs/lock/humann3_env.linux-64.txt"
params.decontam_env = "${projectDir}/conda_envs/lock/decontam_env.linux-64.txt"
params.validation_env = "${projectDir}/conda_envs/lock/validation_env.linux-64.txt"
// Batch correction runs in a PRE-BUILT env (lockfile + cqrReg 1.2.1 + pinned ConQuR, built by
// conda_envs/build_envs.sbatch): a fresh build from the lockfile alone would lack ConQuR (audit D04).
params.batch_corr_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/envs/batch_correction_20260923"
params.antismash_env = "${projectDir}/conda_envs/antismash_env.yml"
params.batch_corr_method = "tune"  // ConQuR mode: "tune" (auto-phased Tune_ConQuR) or "vanilla"


// ============================================================================
// SCRIPT PATHS
// ============================================================================

params.scripts ="${projectDir}/scripts"
params.decontam_script = "${projectDir}/scripts/decontamination.R"
params.batch_corr_script = "${projectDir}/scripts/2500703_batch_correction_normalization.r"

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
include { multiqc } from './Modules/multiqc.nf'
include { qcSummary } from './Modules/qc_summary.nf'
include { cohortReport } from './Modules/cohort_report.nf'
include { mapReads } from './Modules/map_reads.nf'
// Taxonomic classification modules
include { Bracken } from './Modules/Bracken.nf'
include { metaphlan4 } from './Modules/metaphlan4.nf'
include { process_metaphlan; process_bracken; consensus_taxa } from './Modules/preprocess_taxa.nf'
include { humann3; merge_humann3 } from './Modules/humann3.nf'

// Preprocessing modules (optional)
include { Decontamination } from './Modules/decontamination.nf'
include { validationSpeciesList; buildValidationRef; alignForValidation; applyValidationMask; compareDecontamLevels } from './Modules/alignment_validation.nf'
include { BatchCorrection } from './Modules/batch_correction.nf'

// ============================================================================
// PARAMETER VALIDATION HELPERS
// ============================================================================

def pathParam(String name, value) {
    if (value == null) return null
    if (value instanceof Boolean || value.toString().trim() in ['true', 'false']) {
        error "--${name} expects a path but was given as a bare flag (value: '${value}'). Pass --${name} <path>, or omit it."
    }
    def resolved = file(value.toString().trim())
    if (!resolved.exists()) error "--${name} does not exist: ${resolved}"
    return resolved
}

// Dependency digest for a process's `val deps` input (audit C03). Nextflow hashes a task's
// unrendered script and the VALUES of params it references, but not the contents behind a path
// string, so without this an edited helper script or a rebuilt database was served stale on
// -resume. code: content SHA-256 (small files); data: name+size+mtime (databases, indexes).
def depsOf(List code, List data = []) {
    (code.collect { Provenance.code(it) } + data.collect { Provenance.data(it) }).join(';')
}

def requirePathParam(String name, value, String why) {
    if (!value) error "--${name} is required ${why}"
    return pathParam(name, value)
}

// ============================================================================
// BOOLEAN PARAMETER NORMALISATION
// ============================================================================
// Nextflow 24.10 coerced `--flag false` on the command line into a Boolean.
// 26.x does not: the value arrives as the String "false", which is truthy in
// Groovy, so every boolean flag passed as `false` would silently be ON.
// Verified on 26.04.6 (String, truthy) vs 24.10.0 (Boolean, falsy).
// Params cannot be normalised in place: Nextflow ignores every assignment to a
// param after the first ("defined multiple times -- Assignments following the
// first are ignored"). So the entry workflow resolves each flag into a local via
// boolParam() below, and reads the local rather than the param.

def requireRoleParam(String name, value) {
    // A bare flag arrives as Boolean true, the same trap the path params hit.
    if (value == null || value == true || value.toString().trim() == '') {
        error "--${name} is required when decontamination runs. Pass it explicitly, " +
              "or use --cohort_preset tcga for the TCGA sample-type conventions."
    }
    return value.toString()
}

def boolParam(String name, value) {
    if (value instanceof Boolean) return value
    if (value == null) return false
    def s = value.toString().trim()
    // Only true/false, lower case. Modules re-read some flags as `params.x.toString() == 'true'`, so
    // accepting yes/on/1 -- or TRUE -- here made those flags look on at validation but stay OFF in
    // the module (audit C19).
    if (s == 'true')          return true
    if (s in ['false', ''])   return false
    error "--${name} expects a boolean (true/false) but was given '${value}'"
}


// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    // Resolved once: see BOOLEAN PARAMETER NORMALISATION above.
    def skipHostDepletion  = boolParam('skip_host_depletion', params.skip_host_depletion)
    def skipClassification = boolParam('skip_classification', params.skip_classification)
    def skipConsensus      = boolParam('skip_consensus', params.skip_consensus)
    def skipHumann3        = boolParam('skip_humann3', params.skip_humann3)
    def runDecontam        = boolParam('run_decontam', params.run_decontam)
    def runBatchCorrection = boolParam('run_batch_correction', params.run_batch_correction)
    def skipMultiqc        = boolParam('skip_multiqc', params.skip_multiqc)
    def skipCohortReport   = boolParam('skip_cohort_report', params.skip_cohort_report)
    def runValidation      = boolParam('run_alignment_validation', params.run_alignment_validation)
    boolParam('decontam_allow_unknown_library', params.decontam_allow_unknown_library)
    boolParam('batch_corr_allow_uncorrected', params.batch_corr_allow_uncorrected)
    boolParam('save_intermediates', params.save_intermediates)
    boolParam('tumor_only', params.tumor_only)
    boolParam('validation_export_reads', params.validation_export_reads)

    // Channels feeding the cohort QC summary. Collected as the branches run, because
    // which artefacts exist depends on which stages were asked for.
    def qcParts = []


    // ============================================================================
    // PARAMETER VALIDATION
    // ============================================================================
    // Every path parameter is resolved here, before a single task is submitted.
    // A bare `--pangenome_db` (flag with no value) arrives as the Boolean true and
    // used to fail inside mapReads with "Pangenome index path does not exist: true",
    // after host depletion had already run for hours.

    if (!(params.batch_corr_method in ['tune', 'vanilla'])) {
        error "--batch_corr_method must be 'tune' or 'vanilla' (got '${params.batch_corr_method}'); ComBat is not implemented"
    }
    if (!(params.sample_failure_strategy in ['finish', 'ignore'])) {
        error "--sample_failure_strategy must be one of: finish, ignore"
    }
    def adapterTrimMode = params.adapter_trim.toString().toLowerCase()
    if (!(adapterTrimMode in ['auto', 'always', 'never'])) {
        error "Invalid --adapter_trim '${params.adapter_trim}'; expected auto, always or never"
    }

    if (params.start_from == 'beginning' && !skipHostDepletion) {
        requirePathParam('hg38_db', params.hg38_db, 'for host depletion')
        requirePathParam('t2t_phix_db', params.t2t_phix_db, 'for host depletion')
        if (adapterTrimMode != 'never') requirePathParam('adapters', params.adapters, 'for adapter trimming (or use --adapter_trim never)')

        def pangenomeIndex = pathParam('pangenome_db', params.pangenome_db)
        if (pangenomeIndex && pangenomeIndex.isDirectory() && !pangenomeIndex.listFiles().any { it.name.endsWith('.mmi') }) {
            error "--pangenome_db is a directory with no .mmi files in it: ${pangenomeIndex}"
        }
    }

    if (runValidation) {
        // Validation aligns each sample's host-depleted reads, so it needs the full run.
        if (params.start_from != 'beginning' || skipClassification || skipHostDepletion) {
            error "--run_alignment_validation needs --start_from beginning with host depletion and " +
                  "classification on: it aligns each sample's host-depleted reads"
        }
        if (params.validation_refseq_snapshot && !(params.validation_refseq_snapshot.toString() ==~ /\d{4}-\d{2}/)) {
            error "--validation_refseq_snapshot must be YEAR-MONTH, e.g. 2026-09 (got '${params.validation_refseq_snapshot}')"
        }
        if (!params.validation_ref && !params.validation_taxdb && !params.kraken_db) {
            error "--run_alignment_validation needs --kraken_db (for its taxDB) or --validation_taxdb"
        }
        [validation_min_bins_ratio: params.validation_min_bins_ratio,
         validation_min_identity: params.validation_min_identity,
         validation_min_unique_frac: params.validation_min_unique_frac,
         validation_max_missing_frac: params.validation_max_missing_frac].each { name, v ->
            def x = v.toString().isNumber() ? v.toString().toDouble() : -1
            if (x < 0 || x > 1) error "--${name} must be a number in [0, 1] (got '${v}')"
        }
        [validation_min_support: params.validation_min_support,
         validation_min_confirmed: params.validation_min_confirmed,
         validation_max_genomes: params.validation_max_genomes,
         validation_batch_size: params.validation_batch_size].each { name, v ->
            if (!v.toString().isInteger() || v.toString().toInteger() < 1)
                error "--${name} must be a positive integer (got '${v}')"
        }
        if (params.validation_groups && !file(params.validation_groups).exists())
            error "--validation_groups does not exist: ${params.validation_groups}"
    }

    // Integer-only: the Bracken gate compares with bash [[ -ge ]], and a value like 1e5 made that
    // test error out so every sample silently skipped MetaPhlAn (audit A13).
    ['consensus_min_reads': params.consensus_min_reads, 'bracken_threshold': params.bracken_threshold].each { name, v ->
        if (!v.toString().isInteger() || v.toString().toInteger() < 0)
            error "--${name} must be a non-negative integer (got '${v}')"
    }

    if (params.start_from == 'beginning' && !skipClassification) {
        requirePathParam('kraken_db', params.kraken_db, 'for KrakenUniq/Bracken classification')
        requirePathParam('metaphlan_db', params.metaphlan_db, 'for MetaPhlAn 4 classification')
        // Contents, not just existence (audit C17): the files these steps actually open.
        def krakenDir = file(params.kraken_db)
        if (krakenDir.isDirectory()) {
            ['taxDB', "database${params.bracken_read_length}mers.kmer_distrib"].each { f ->
                if (!file("${krakenDir}/${f}").exists())
                    error "--kraken_db ${krakenDir} has no ${f}" +
                          (f.endsWith('kmer_distrib') ? " (Bracken needs one built for --bracken_read_length ${params.bracken_read_length})" : "")
            }
        }
        def mpaDir = file(params.metaphlan_db)
        if (mpaDir.isDirectory() && !file("${mpaDir}/${params.metaphlan_index}.pkl").exists()) {
            error "--metaphlan_db ${mpaDir} has no ${params.metaphlan_index}.pkl (check --metaphlan_index)"
        }
        if (!skipHumann3) {
            requirePathParam('humann3_nucleotide_db', params.humann3_nucleotide_db, 'for HUMAnN3 functional profiling')
            requirePathParam('humann3_protein_db', params.humann3_protein_db, 'for HUMAnN3 functional profiling')
        }
    }

    if (!(params.start_from in ['beginning', 'decontam', 'batch_correction'])) {
        error "--start_from must be one of: beginning, decontam, batch_correction"
    }

    // --- Decontamination sample-role contract -------------------------------------
    // Scoped to runs that actually reach decontamination. June-young's tree validates
    // --taxonomic_rank unconditionally, which would reject every host-depletion-only
    // run; that is the one change made while porting his checks.
    // Only when decontam actually runs: --start_from decontam --run_decontam false must reach
    // the named C20 error ("Batch correction needs an input table"), not a rank error for a
    // step that is switched off.
    def needsDecontamContract = runDecontam && params.start_from != 'batch_correction'

    if (params.cohort_preset && params.cohort_preset != 'tcga') {
        error "Invalid --cohort_preset '${params.cohort_preset}'; expected tcga"
    }
    if (params.cohort_preset == 'tcga') {
        if (params.control_role != 'surrogate_biological') {
            error "--cohort_preset tcga conflicts with --control_role '${params.control_role}'"
        }
        if (params.positive_values != 'Primary Tumor,Solid Tissue Normal') {
            error "--cohort_preset tcga conflicts with --positive_values '${params.positive_values}'"
        }
        if (params.control_values != 'Blood Derived Normal') {
            error "--cohort_preset tcga conflicts with --control_values '${params.control_values}'"
        }
        if (params.type_column != 'Sample Type_x') {
            error "--cohort_preset tcga conflicts with --type_column '${params.type_column}'"
        }
    }

    if (needsDecontamContract) {
        if (!params.taxonomic_rank) {
            error "--taxonomic_rank is required when decontamination runs; expected genus or species"
        }
        if (!(params.taxonomic_rank in ['genus', 'species'])) {
            error "--taxonomic_rank must be one of: genus, species (got '${params.taxonomic_rank}')"
        }
        def controlRole = requireRoleParam('control_role', params.control_role)
        if (!(controlRole in ['experimental_blank', 'surrogate_biological'])) {
            error "Invalid --control_role '${controlRole}'; expected experimental_blank or surrogate_biological"
        }
        requireRoleParam('positive_values', params.positive_values)
        requireRoleParam('control_values', params.control_values)
        // The metadata must carry the columns decontam will use (audit C17): fail at launch, not
        // after classification has run.
        if (params.metadata_file && file(params.metadata_file).exists()) {
            def header = file(params.metadata_file).withReader { it.readLine() } ?: ''
            def cols = header.split(header.contains('\t') ? '\t' : ',').collect { it.trim() }
            [params.sample_id_column, params.type_column, params.batch_column].findAll { it }.each { c ->
                if (!(c in cols)) error "--metadata_file has no '${c}' column (columns: ${cols.join(', ')})"
            }
        }
        if (!(params.unselected_type_policy in ['error', 'exclude'])) {
            error "--unselected_type_policy must be one of: error, exclude"
        }
        if (!params.decontam_min_library.toString().isNumber() || params.decontam_min_library.toString().toDouble() < 0) {
            error "--decontam_min_library must be a number >= 0"
        }
        if (!(params.decontam_unevaluable_policy in ['skip', 'error'])) {
            error "--decontam_unevaluable_policy must be one of: skip, error"
        }
        if (!(params.decontam_zero_total_policy in ['keep', 'drop'])) {
            error "--decontam_zero_total_policy must be one of: keep, drop"
        }
        // Pre-existing hole: with classification and consensus both skipped there is no
        // Bracken channel to feed decontamination, and the workflow failed at runtime with
        // "No such variable: BRACKEN_FILES". Say so up front instead.
        if (params.start_from == 'beginning' && skipConsensus && skipClassification) {
            error "Cannot run decontamination with both --skip_classification and " +
                  "--skip_consensus from the beginning: there is no taxonomic table to " +
                  "decontaminate. Either enable classification, or start downstream with " +
                  "--start_from decontam --consensus_otu_table <path>."
        }
    }

    // Consensus consumes classification output; without it the workflow died with an unnamed
    // "No such variable" error (Workflow 2 critic).
    if (params.start_from == 'beginning' && skipClassification && !skipConsensus) {
        error "--skip_consensus true is required with --skip_classification: consensus consumes the classification output"
    }

    // Batch correction's own contract, checked at launch rather than when the task starts after
    // every upstream stage (audit C04/C17; decision 16 for the covariates).
    if (runBatchCorrection) {
        // No input table for batch correction: decided HERE, before any channel exists. Raised
        // later inside the workflow, the error raced the sample-sheet reader, which could surface
        // a misleading "Sample sheet has no data rows" first (audit C20).
        def hasBatchInput = params.start_from == 'batch_correction' || runDecontam ||
                            (params.start_from == 'beginning' && (!skipConsensus || !skipClassification))
        if (!hasBatchInput) {
            error "Batch correction needs an input table: run it after decontamination (--run_decontam true), " +
                  "after consensus (--start_from beginning), or use --start_from batch_correction with " +
                  "--decontam_otu_table and --metadata_file (audit C20)"
        }
        if (!params.taxonomic_rank || !(params.taxonomic_rank in ['genus', 'species'])) {
            error "--taxonomic_rank must be genus or species for batch correction (got '${params.taxonomic_rank ?: ''}')"
        }
        // An empty or short pin would make the module's startsWith() check accept any install (audit D04).
        if (params.conqur_sha.toString().trim().length() < 7) {
            error "--conqur_sha must be a commit SHA of at least 7 characters (got '${params.conqur_sha}')"
        }
        if (params.phase.toString() != 'auto') {
            error "--phase: only 'auto' is implemented (got '${params.phase}')"
        }
        def covariates = params.batch_corr_covariates.toString().split(',').collect { it.trim() }.findAll { it }
        if (!covariates || covariates.collect { it.toLowerCase() } == ['none']) {
            error "ConQuR needs >= 1 covariate that varies within batches; --batch_corr_covariates " +
                  "'${params.batch_corr_covariates}' cannot be batch-corrected (decision 16). " +
                  "Name the covariates, or pass --run_batch_correction false."
        }
        if (!params.metadata_file) {
            error "--metadata_file is required for batch correction"
        }
        def metaFile = file(params.metadata_file)
        if (!metaFile.exists()) error "--metadata_file does not exist: ${params.metadata_file}"
        def header = metaFile.withReader { it.readLine() } ?: ''
        def cols = header.split(header.contains('\t') ? '\t' : ',').collect { it.trim() }
        def needed = [params.sample_id_column, params.batch_column] + covariates
        if (params.tumor_only.toString() == 'true') needed << params.type_column
        def missing = needed.findAll { it && !(it in cols) }.unique()
        if (missing) {
            error "Batch correction needs metadata column(s) that --metadata_file lacks: ${missing.join(', ')} " +
                  "(columns: ${cols.join(', ')}). Fix --batch_corr_covariates / --batch_column / " +
                  "--sample_id_column, or pass --run_batch_correction false."
        }
    }
    // A flag combination that selects no stage used to build an empty DAG and report Success
    // (Workflow 2 critic): say so instead.
    def anyStage = (params.start_from == 'beginning' && (!skipHostDepletion || !skipClassification || !skipConsensus)) ||
                   runDecontam || runBatchCorrection
    if (!anyStage) {
        error "Nothing to run: --start_from ${params.start_from} with these flags selects no stage " +
              "(host depletion, classification, consensus, decontamination and batch correction are all off)"
    }
    if (params.start_from == 'decontam' && params.consensus_otu_table &&
        !file(params.consensus_otu_table).exists()) {
        error "--consensus_otu_table does not exist: ${params.consensus_otu_table}"
    }
    if (params.start_from == 'batch_correction' && params.decontam_otu_table &&
        !file(params.decontam_otu_table).exists()) {
        error "--decontam_otu_table does not exist: ${params.decontam_otu_table}"
    }

    // ============================================================================
    // PROVENANCE -- launch half
    // ============================================================================
    // Written before any task is submitted, so a run that dies still records what it
    // was attempting. scripts/finalize_provenance.py fills in the outcome, the conda
    // versions actually installed and the output digests, from workflow.onComplete in
    // main.nf. Ported from pksProfiler F16.
    //
    // Reference databases are fingerprinted by name/size/mtime (Provenance.data), not
    // by content: these are multi-gigabyte indexes and hashing them would dominate the
    // run. Scripts are fingerprinted by content (Provenance.code), because that is the
    // claim the record is actually making. The sample sheet and metadata are small and
    // are content-hashed (Provenance.content); only the SHA-256 enters the record (audit C06).
    def runFlags = [
        skip_host_depletion : skipHostDepletion,
        skip_classification : skipClassification,
        skip_consensus      : skipConsensus,
        skip_humann3        : skipHumann3,
        skip_multiqc        : skipMultiqc,
        run_decontam        : runDecontam,
        run_batch_correction: runBatchCorrection,
        run_alignment_validation: runValidation,
        skip_cohort_report  : skipCohortReport,
    ]
    // The taxonomy DB alignment validation actually reads: mirrors STEP 3b (audit C06).
    def validationTaxdb = params.validation_ref ? null :
        (params.validation_taxdb ?: (params.kraken_db ? "${params.kraken_db}/taxDB" : null))
    RunRecord.write(
        // One record per LAUNCH: `-resume` reuses the session id, so keying on it alone overwrote
        // the previous launch's record (audit C06). The run name is new for every launch.
        RunRecord.path(params.tracedir, "${workflow.sessionId}_${workflow.runName}"),
        RunRecord.start([
            session_id  : workflow.sessionId.toString(),
            code        : RunRecord.gitIdentity(workflow.projectDir, workflow.commitId, workflow.revision) + [
                repository      : workflow.manifest.name,
                project_dir     : workflow.projectDir.toString(),
                script_id       : workflow.scriptId?.toString(),
                scripts_digest  : Provenance.code(params.scripts),
                modules_digest  : Provenance.code("${projectDir}/Modules"),
                decontam_script : Provenance.code(params.decontam_script),
                batch_corr_script: Provenance.code(params.batch_corr_script),
                envs_digest     : Provenance.specs("${projectDir}/conda_envs"),   // yml + lock/ (audit D08)
            ],
            invocation  : [
                command_line    : workflow.commandLine,
                launch_dir      : workflow.launchDir.toString(),
                work_dir        : workflow.workDir.toString(),
                profile         : workflow.profile,
                sample_sheet    : [path: params.sample?.toString(), digest: Provenance.content(params.sample)],
                metadata_file   : [path: params.metadata_file?.toString(), digest: Provenance.content(params.metadata_file)],
                stages_requested: RunRecord.stages(params, runFlags),
                params          : params.collectEntries { key, value -> [(key): value?.toString()] },
            ],
            dependencies: [
                hg38_db              : Provenance.data(params.hg38_db),
                t2t_phix_db          : Provenance.data(params.t2t_phix_db),
                pangenome_db         : Provenance.data(params.pangenome_db),
                kraken_db            : Provenance.data(params.kraken_db),
                metaphlan_db         : Provenance.data(params.metaphlan_db),
                humann3_nucleotide_db: Provenance.data(params.humann3_nucleotide_db),
                humann3_protein_db   : Provenance.data(params.humann3_protein_db),
                cram_reference       : Provenance.data(params.cram_reference),
                adapters             : Provenance.code(params.adapters),
                // Alignment validation lane (audit C06). Recorded only when it runs: the cache
                // listing is not free on Lustre, and an unused path is not a dependency.
                validation_ref       : runValidation ? Provenance.data(params.validation_ref) : 'unused',
                validation_taxdb     : runValidation ? Provenance.data(validationTaxdb) : 'unused',
                validation_groups    : runValidation ? Provenance.code(params.validation_groups) : 'unused',
                validation_cache_dir : runValidation ? Provenance.data(params.validation_cache_dir) : 'unused',
            ],
            environment : [
                nextflow_version: workflow.nextflow.version?.toString(),
                nextflow_build  : workflow.nextflow.build?.toString(),
                user_name       : workflow.userName?.toString(),
                conda_enabled   : workflow.containerEngine == null,
                conda_cache_dir : (System.getenv('NXF_CONDA_CACHEDIR') ?: "${projectDir}/.conda_cache"),
                // Spec of each env the requested stages run in; the finaliser digests these
                // rather than walking conda_cache_dir (audit C06/D08).
                conda_envs      : RunRecord.condaEnvs(params, runFlags),
            ],
        ])
    )
    log.info "Provenance record: ${params.tracedir}/provenance/${workflow.sessionId}_${workflow.runName}.json"

    // ============================================================================
    // CONDITIONAL WORKFLOW BASED ON ENTRY POINT
    // ============================================================================

    if (params.start_from == 'beginning') {
        // Full pipeline from BAM files

        if (skipHostDepletion && !skipClassification) {
            error "--skip_host_depletion requires --skip_classification true; classification consumes mapped reads produced by Module 1"
        }

        // ------------------- STEP1: HOST DEPLETION ---------------------- //

        if (!skipHostDepletion) {
            // Read and parse the sample sheet
            def requestedInputType = params.input_data_type.toString().toLowerCase()
            if (!(requestedInputType in ['auto', 'bam', 'cram', 'fastq'])) {
                error "Invalid --input_data_type '${params.input_data_type}'; expected auto, bam, cram, or fastq"
            }

            sample_sheet = nextflow.Channel.fromPath(params.sample, checkIfExists: true)
                .splitCsv(header: true)
                .map { row ->
                    def sampleID = row.patient?.trim()
                    if (!sampleID) error "Sample sheet row is missing a non-empty 'patient' value"
                    // IDs become file names, awk variables, CLI lists and table columns downstream.
                    // Characters outside this set broke the validation batch (B07) and merges (C16).
                    if (!(sampleID ==~ /[A-Za-z0-9._-]+/)) {
                        error "Sample ID '${sampleID}' contains characters outside [A-Za-z0-9._-]; rename it in the sample sheet"
                    }
                    def bam = row.bam?.trim()
                    def cram = row.cram?.trim()
                    def r1 = (row.fastq_1 ?: row.r1)?.trim()
                    def r2 = (row.fastq_2 ?: row.r2)?.trim()
                    def presentTypes = []
                    if (bam) presentTypes << 'bam'
                    if (cram) presentTypes << 'cram'
                    if (r1 || r2) presentTypes << 'fastq'
                    if (requestedInputType == 'auto' && presentTypes.size() > 1) {
                        error "Sample '${sampleID}' has ambiguous inputs (${presentTypes.join(', ')}); populate only one input type"
                    }
                    def inputType = requestedInputType == 'auto' ? (presentTypes ? presentTypes[0] : null) : requestedInputType
                    if (!inputType) error "Sample '${sampleID}' has no input. Provide bam, cram, or fastq_1/r1"
                    if (inputType == 'bam' && !bam) error "Sample '${sampleID}' requires a value in the bam column"
                    if (inputType == 'cram' && !cram) error "Sample '${sampleID}' requires a value in the cram column"
                    if (inputType == 'fastq' && !r1) error "Sample '${sampleID}' requires fastq_1 (or r1); fastq_2/r2 is optional"
                    if (inputType == 'fastq' && r2 && r1 == r2) error "Sample '${sampleID}' has the same path for both FASTQ mates"
                    def reference = (row.cram_reference ?: params.cram_reference)?.toString()?.trim()
                    def inputs = inputType == 'bam' ? [file(bam, checkIfExists: true)] :
                        (inputType == 'cram' ? [file(cram, checkIfExists: true)] :
                        ([file(r1, checkIfExists: true)] + (r2 ? [file(r2, checkIfExists: true)] : [])))
                    def references = inputType == 'cram' && reference ? [file(reference, checkIfExists: true)] : []
                    tuple(sampleID, inputType, inputs, references)
                }

            // A checkpoint every row must pass through, so the check always runs.
            // Two rows with the same patient ID -- or IDs differing only by surrounding
            // whitespace, which `trim()` collapses -- used to flow straight through and
            // then collide in publishDir and in every merged table, one sample silently
            // overwriting the other. This is the same defect as pksProfiler F03.
            sample_sheet = sample_sheet
                .toList()
                .map { rows ->
                    if (!rows) {
                        error "Sample sheet has no data rows: ${params.sample}"
                    }
                    def ids = rows.collect { it[0] }
                    def duplicates = ids.countBy { it }.findAll { id, n -> n > 1 }.keySet()
                    if (duplicates) {
                        error "Duplicate sample IDs in the sample sheet: " +
                              "${duplicates.sort().join(', ')}. Sample IDs must be unique " +
                              "-- note that surrounding whitespace is stripped, so ' A ' " +
                              "and 'A' are the same ID."
                    }
                    return rows
                }
                .flatMap()

            // Every sample the sheet asked for, so the QC summary can name one that produced
            // nothing -- e.g. dropped by --sample_failure_strategy ignore (audit C10).
            qcParts << sample_sheet.map { it[0] }.collectFile(name: 'expected_samples.txt', newLine: true, sort: true)

            sample_sheet.branch {
                alignment: it[1] in ['bam', 'cram']
                fastq: it[1] == 'fastq'
            }.set { INPUTS_BY_TYPE }

            extractReads(INPUTS_BY_TYPE.alignment.map { sampleID, inputType, inputs, references -> tuple(sampleID, inputType, inputs[0], references) },
                         depsOf(["${params.scripts}/validate_cram_reference.py"], [params.cram_reference]))
            extractReads.out.reads.set { EXTRACTED_READS }
            qcParts << extractReads.out.counts   // *.extract_counts.tsv: BAM/CRAM library denominator (C07)
            // `pretrimmed` marks reads that came out of an aligned BAM/CRAM: whoever
            // produced the alignment adapter-trimmed them before aligning.
            INPUTS_BY_TYPE.fastq
                .map { sampleID, inputType, reads, references -> tuple(sampleID, reads, false) }
                .mix(EXTRACTED_READS.map { sampleID, reads -> tuple(sampleID, reads, true) })
                .set { UNMAPPED_READS }

            // Perform FastQC on the extracted fastq files
            FASTQC1(UNMAPPED_READS.map { sampleID, reads, pretrimmed -> tuple(sampleID, reads) })

            // adapter_trim: `auto` scans raw FASTQ against --adapters and skips that scan for
            // reads extracted from an alignment. Measured on an aligned CRAM: fastp with
            // --adapter_fasta (234 sequences) took 2 h 50 m and trimmed 0 reads, against
            // 3 m 15 s without it. `always` and `never` override the per-sample decision.
            // Dropping --adapter_fasta does not disable adapter trimming: fastp still
            // auto-detects adapters.
            UNMAPPED_READS
                .map { sampleID, reads, pretrimmed ->
                    def trimAdapters = adapterTrimMode == 'always' ? true :
                        (adapterTrimMode == 'never' ? false : !pretrimmed)
                    tuple(sampleID, reads, trimAdapters)
                }
                .set { FILTER_INPUT }

            // Filter poor quality reads using fastp
            filterReads(FILTER_INPUT, depsOf([params.adapters]))
            filterReads.out.reads.set { FILTERED_ALL }

            // A sample that has no reads left after filtering is dropped with a warning
            // instead of failing the run. Aligned inputs with no unmapped records (common
            // for RNA-seq BAMs written without --outSAMunmapped Within) land here.
            FILTERED_ALL
                .branch { sampleID, reads, readCount ->
                    pass: readCount.text.trim().toLong() > 0
                    empty: true
                }
                .set { FILTERED_BY_DEPTH }

            FILTERED_BY_DEPTH.empty.subscribe { sampleID, reads, readCount ->
                log.warn "Sample ${sampleID} has 0 reads after filtering; dropping it from the rest of the run"
            }

            FILTERED_BY_DEPTH.pass
                .map { sampleID, reads, readCount -> tuple(sampleID, reads) }
                .set { FILTERED_UNMAPPED_READS }

            // Perform FastQC on the filtered fastq files
            FASTQC2(FILTERED_UNMAPPED_READS)

            mapReads(FILTERED_UNMAPPED_READS, depsOf([], [params.hg38_db, params.t2t_phix_db, params.pangenome_db]))
            mapReads.out.reads.set { MAPPED_READS }

            MAPPED_READS.multiMap { sampleID, hg38Reads, t2tReads, hostDepletedReads ->
                Hg38: tuple(sampleID, hg38Reads)
                T2T: tuple(sampleID, t2tReads)
                PAN: tuple(sampleID, hostDepletedReads)
            }
            .set { MAPPED_READS_MULTI }

            // Perform FastQC on the mapped Reads
            FASTQCHG38(MAPPED_READS_MULTI.Hg38)
            FASTQCT2T(MAPPED_READS_MULTI.T2T)
            FASTQCPANGENOME(MAPPED_READS_MULTI.PAN)

            // Cohort QC report. fastp's per-sample JSON carries reads in/out, duplication
            // and adapter stats; the five FastQC passes cover raw, filtered, and each
            // depletion stage. `.collect()` makes this one task over the whole cohort.
            if (!skipMultiqc) {
                Channel.empty()
                    .mix(filterReads.out.fastp_json)
                    .mix(FASTQC1.out.zip)
                    .mix(FASTQC2.out.zip)
                    .mix(FASTQCHG38.out.zip)
                    .mix(FASTQCT2T.out.zip)
                    .mix(FASTQCPANGENOME.out.zip)
                    .collect()
                    .set { QC_FILES }
                multiqc(QC_FILES)
            } else {
                log.info "Skipping MultiQC cohort report"
            }

            // Feed the cohort QC summary regardless of --skip_multiqc: the summary table
            // and the MultiQC report answer different questions.
            qcParts << filterReads.out.fastp_json
            qcParts << mapReads.out.counts
        } else {
            log.info "Skipping host depletion step"
        }

        // ------------------- STEP2: TAXONOMIC CLASSIFICATION ---------------------- //

        if (!skipClassification) {
            // BRACKEN
            Bracken(MAPPED_READS_MULTI.PAN, depsOf(["${params.scripts}/kreport2mpa.py"], [params.kraken_db]))
            Bracken.out.reports.set { BRACKEN_OUT }

            // Collect all kraken reports once all samples are done
            BRACKEN_OUT.map { sampleID, kraken_report, classified_fasta, unclassified_fasta, bracken_reports, bracken_krakenreports, bracken_mpa_reports ->
                tuple(bracken_mpa_reports)
            }
            .flatten()
            .collect()
            .set {Bracken_mpa_files}

            // Run process_bracken once all reports are available
            process_bracken(Bracken_mpa_files, depsOf(["${params.scripts}/combine_mpa.py"])).set { BRACKEN_FILES }

            // QC summary inputs from classification: the KrakenUniq report carries the
            // microbial read count the consensus gate uses, and the Bracken genus/species
            // reports give the per-sample taxon counts.
            qcParts << BRACKEN_OUT.map { it[1] }          // *.krakenuniq.report.txt
            qcParts << BRACKEN_OUT.flatMap { it[4] }      // *.bracken.{G,S}.report.txt

            // METAPHLAN4: only samples passing the per-sample Bracken read gate
            // Only samples that cleared --consensus_min_reads emit this channel at all,
            // so no filtering is needed: absence is the gate.
            Bracken.out.metaphlan_reads.set { METAPHLAN_INPUT }
            metaphlan4(METAPHLAN_INPUT, depsOf([], [params.metaphlan_db])).set { METAPHLAN_OUT }

            // Collect all MetaPhlAn profiles once all samples are done
            METAPHLAN_OUT.map { sampleID, bowtie2_file, sam_file, profiled_metagenome ->
                tuple(profiled_metagenome)
            }
            .flatten()
            .collect()
            .set {metaphlan4_files}

            METAPHLAN_OUT
                .map { sampleID, bowtie2_file, sam_file, profiled_metagenome -> tuple(sampleID, profiled_metagenome) }
                .set { METAPHLAN_PROFILES }

            process_metaphlan(metaphlan4_files, depsOf(["${params.scripts}/merge_metaphlan_tables.py"])).set { METAPHLAN_FILES }

            qcParts << METAPHLAN_OUT.map { it[3] }        // *.profiled_metagenome.txt

            // HUMANN3 (functional profiling - runs in parallel with classification merging)
            if (!skipHumann3) {
                // HUMAnN3 runs only on samples MetaPhlAn ran on, i.e. those that cleared
                // --consensus_min_reads, and uses that sample's own profile. Below the gate,
                // HUMAnN would find almost no species for its nucleotide tier and fall back
                // to a translated search of every read: its slowest mode, and the one most
                // exposed to residual human reads. Functional profiling needs more depth
                // than taxonomy, not less. Decided 2026-09-23 (Ammal). Which samples ran is
                // visible per sample in the QC summary (metaphlan_genera is NA when skipped).
                // Taking the profile from the channel rather than the published METAPHLAN4
                // directory avoids racing metaphlan4 itself.
                MAPPED_READS_MULTI.PAN
                    .join(METAPHLAN_PROFILES)
                    .set { HUMANN3_INPUT }
                HUMANN3_INPUT.count().subscribe { n ->
                    if (n == 0) log.warn "HUMAnN3 was requested but no sample cleared " +
                        "--consensus_min_reads (${params.consensus_min_reads}); no functional profiles will be produced."
                }

                humann3(HUMANN3_INPUT, depsOf([], [params.humann3_nucleotide_db, params.humann3_protein_db, params.metaphlan_db])).set { HUMANN3_OUT }

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
        if (!skipClassification) {
            BRACKEN_FILES
                .map { bracken_genus_file, bracken_species_file -> tuple(bracken_genus_file, bracken_species_file) }
                .set { bracken_file_tuple }
        }

        if (!skipConsensus) {
            // Step 1: Get merged table
            // If no sample cleared --consensus_min_reads, metaphlan4 never runs and the
            // collect() above emits nothing. Without a stand-in, consensus_taxa would then
            // silently never run and the workflow would report Success with no table for
            // decontam. The sentinel becomes [] below, and the script passes Bracken through.
            METAPHLAN_FILES
                .map { merged_table, merged_genus, merged_species, merged_SGB -> merged_genus }
                .ifEmpty('NO_METAPHLAN')
                .set { metaphlan4_merged_table }

            // Step 2: Combine for consensus (bracken_file_tuple already defined above)
            metaphlan4_merged_table
                .combine(bracken_file_tuple)
                .map { mp, genus, species -> tuple(mp == 'NO_METAPHLAN' ? [] : mp, genus, species) }
                .set { consensus_input }

            // consensus_taxa emits three outputs (pdf, common genus, common species).
            // `.set` on the multi-channel object made every downstream `.map` fail with
            // "Multi-channel output cannot be applied to operator map". Select the one
            // table the requested rank needs, and take it from the channel rather than
            // re-reading it out of publishDir (which raced with the publish).
            consensus_taxa(consensus_input, depsOf(["${params.scripts}/compute_consensus_taxa.py"]))
            if (params.taxonomic_rank == 'species') {
                consensus_taxa.out[2].set { CONSENSUS_OUTPUT }
            } else {
                consensus_taxa.out[1].set { CONSENSUS_OUTPUT }
            }
            consensus_taxa.out[4].set { CONSENSUS_LIBRARY_SIZES }
            consensus_taxa.out[3].set { CONSENSUS_STATUS }
        } else {
            log.info "Skipping consensus taxa step"
        }

        // -------------------  STEP 3b: ALIGNMENT VALIDATION (OPTIONAL) ---------------------- //
        if (runValidation) {
            // The species table going downstream: consensus (or its Bracken pass-through), else Bracken.
            // `.set` rather than `def x = <channel>`: Nextflow 24.10 rejects the latter
            // ("Variable already defined in the process scope").
            if (!skipConsensus) {
                consensus_taxa.out[2].set { VALIDATION_TABLE }
            } else {
                BRACKEN_FILES.map { genus_file, species_file -> species_file }.set { VALIDATION_TABLE }
            }
            def groupsFile = params.validation_groups ? file(params.validation_groups, checkIfExists: true) : []
            def taxdb = params.validation_ref ? [] :
                file(params.validation_taxdb ?: "${params.kraken_db}/taxDB", checkIfExists: true)

            if (params.validation_ref) {
                def refDir = file(params.validation_ref, checkIfExists: true)
                Channel.of(tuple('prebuilt',
                        file("${refDir}/validation_reference.mmi", checkIfExists: true),
                        file("${refDir}/contig_labels.tsv", checkIfExists: true),
                        file("${refDir}/validation_reference_manifest.tsv", checkIfExists: true),
                        file("${refDir}/missing_species.tsv", checkIfExists: true)))
                    .first().set { VALIDATION_REF }
            } else {
                validationSpeciesList(VALIDATION_TABLE, taxdb, depsOf(["${params.scripts}/validation_species_list.py", "${params.scripts}/validation_common.py", "${params.scripts}/build_validation_ref.py"]))
                buildValidationRef(validationSpeciesList.out.species)
                buildValidationRef.out.ref.first().set { VALIDATION_REF }
            }
            MAPPED_READS_MULTI.PAN
                .collate(params.validation_batch_size.toString().toInteger())
                .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }
                .set { VALIDATION_BATCHES }
            alignForValidation(VALIDATION_BATCHES, VALIDATION_REF,
                               VALIDATION_TABLE.first(), groupsFile,
                               depsOf(["${params.scripts}/validation_score.py", "${params.scripts}/validation_common.py"]))
            applyValidationMask(VALIDATION_TABLE, alignForValidation.out.stats.flatten().collect(), groupsFile,
                                depsOf(["${params.scripts}/apply_validation_mask.py", "${params.scripts}/validation_common.py"]))
            applyValidationMask.out.species.set { VALIDATED_SPECIES }
            applyValidationMask.out.genus.set { VALIDATED_GENUS }
            applyValidationMask.out.report.set { VALIDATION_MASK_REPORT }
        }

    } // End of: if (params.start_from == 'beginning')

    // ============================================================================
    // STEP 4: DECONTAMINATION (OPTIONAL)
    // ============================================================================

    if (runDecontam && params.start_from != 'batch_correction') {
        // Check if starting from decontamination step
        if (runValidation) {
            // Decontam runs twice on the validated (masked) tables: species is primary, the genus
            // roll-up is the continuity check. Decided 2026-09-23 (Ammal).
            if (!params.metadata_file) {
                error "ERROR: When using --run_decontam, you must provide:\n  --metadata_file <path>"
            }
            def meta = file(params.metadata_file, checkIfExists: true)
            // Library sizes: consensus's bracken_library_sizes.tsv, or -- consensus skipped -- the
            // UNMASKED Bracken species table itself (its column sums are the libraries). Never an
            // empty [] through .combine(), which flattened it away and crashed the map (audit A10).
            if (!skipConsensus) {
                CONSENSUS_LIBRARY_SIZES.set { VALIDATION_LIB_SIZES }
            } else {
                VALIDATION_TABLE.set { VALIDATION_LIB_SIZES }
            }
            VALIDATED_SPECIES.combine(VALIDATION_LIB_SIZES)
                .map { t, lib -> tuple('validated_species', t, meta, lib, 'species') }
                .mix(VALIDATED_GENUS.combine(VALIDATION_LIB_SIZES)
                     .map { t, lib -> tuple('validated_genus', t, meta, lib, 'genus') })
                .set { ch_for_decontam }
        } else if (params.start_from == 'decontam') {
            // Load data from parameters for decontamination entry point
            if (!params.consensus_otu_table || !params.metadata_file) {
                error "ERROR: When starting from 'decontam', you must provide:\n" +
                      "  --consensus_otu_table <path>\n" +
                      "  --metadata_file <path>"
            }
            Channel.of(tuple(
                'decontam_run',
                file(params.consensus_otu_table),
                file(params.metadata_file),
                params.decontam_library_sizes ? file(params.decontam_library_sizes, checkIfExists: true) : [],
                params.taxonomic_rank
            )).set { ch_for_decontam }
        } else if (params.start_from == 'beginning' && !skipConsensus) {
            if (!params.metadata_file) {
                error "ERROR: When using --run_decontam, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            CONSENSUS_OUTPUT
                .combine(CONSENSUS_LIBRARY_SIZES)
                .map { consensus_file, library_sizes ->
                    tuple('consensus_run', consensus_file,
                          file(params.metadata_file, checkIfExists: true), library_sizes,
                          params.taxonomic_rank)
                }
                .set { ch_for_decontam }
        } else if (params.start_from == 'beginning' && skipConsensus) {
            if (!params.metadata_file) {
                error "ERROR: When using --run_decontam with --skip_consensus, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            log.info "Skipping consensus taxa - using Bracken ${params.taxonomic_rank ?: 'genus'} output directly for decontamination"
            BRACKEN_FILES
                .map { bracken_genus_file, bracken_species_file ->
                    def selected = params.taxonomic_rank == 'species' ? bracken_species_file : bracken_genus_file
                    // the unmasked Bracken table is its own library source (audit A07/A10)
                    tuple('bracken_run', selected, file(params.metadata_file, checkIfExists: true), selected,
                          params.taxonomic_rank)
                }
                .set { ch_for_decontam }
        } else {
            // A named contract error instead of logging and then dereferencing an undefined channel (C20).
            error "Decontamination cannot run with --start_from ${params.start_from}: use --start_from beginning " +
                  "(classification on) or --start_from decontam with --consensus_otu_table and --metadata_file"
        }

        Decontamination(ch_for_decontam, depsOf([params.decontam_script]))
        // Decontam roles/status per sample into the QC summary (audit A09/C14). With validation on,
        // decontam runs at both levels; only the primary level (--taxonomic_rank) reports, so the QC
        // columns are deterministic.
        qcParts << Decontamination.out.decontaminated
            .filter { !runValidation || it[0] == "validated_${params.taxonomic_rank}".toString() }
            .map { it[3] }   // <prefix>.sample_validation.tsv
        if (runValidation) {
            // Batch correction continues from the level --taxonomic_rank asks for.
            Decontamination.out.for_batch_correction
                .filter { it[0] == "validated_${params.taxonomic_rank}" }
                .set { ch_for_batch_corr }
            Decontamination.out.decontaminated
                .map { it -> tuple(it[0], it[6]) }            // (prefix, contaminants_final.tsv)
                .branch {
                    species: it[0] == "validated_species"
                    genus: it[0] == "validated_genus"
                }
                .set { FINALS }
            compareDecontamLevels(FINALS.species.map { it[1] }, FINALS.genus.map { it[1] },
                                  depsOf(["${params.scripts}/compare_decontam_levels.py", "${params.scripts}/validation_common.py"]))
        } else {
            Decontamination.out.for_batch_correction.set { ch_for_batch_corr }
        }
    }

    // ============================================================================
    // STEP 5: BATCH CORRECTION (OPTIONAL)
    // ============================================================================

    if (runBatchCorrection) {
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
        } else if (runDecontam) {
            // Use output from decontamination
            ch_for_batch_corr.set { ch_batch_input }
        } else if (params.start_from == 'beginning' && runValidation) {
            // Validation on, decontam off: correct the MASKED table, not the unmasked consensus
            // (Workflow 2 critic). Validation needs classification + consensus, so both exist.
            if (!params.metadata_file) {
                error "ERROR: When using --run_batch_correction without decontam, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            (params.taxonomic_rank == 'species' ? VALIDATED_SPECIES : VALIDATED_GENUS)
                .map { t -> tuple("validated_${params.taxonomic_rank}".toString(), t,
                                  file(params.metadata_file, checkIfExists: true)) }
                .set { ch_batch_input }
        } else if (params.start_from == 'beginning' && !skipConsensus) {
            // No decontam, use consensus output directly
            if (!params.metadata_file) {
                error "ERROR: When using --run_batch_correction without decontam, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            CONSENSUS_OUTPUT
                .map { consensus_file ->
                    tuple('consensus_run', consensus_file,
                          file(params.metadata_file, checkIfExists: true))
                }
                .set { ch_batch_input }
        } else if (params.start_from == 'beginning' && skipConsensus && !skipClassification) {
            if (!params.metadata_file) {
                error "ERROR: When using --run_batch_correction with --skip_consensus, you must provide:\n" +
                      "  --metadata_file <path>"
            }
            log.info "Skipping consensus taxa - using Bracken ${params.taxonomic_rank ?: 'genus'} output directly for batch correction"
            BRACKEN_FILES
                .map { bracken_genus_file, bracken_species_file ->
                    def selected = params.taxonomic_rank == 'species' ? bracken_species_file : bracken_genus_file
                    tuple('bracken_run', selected, file(params.metadata_file, checkIfExists: true))
                }
                .set { ch_batch_input }
        } else {
            error "Batch correction needs an input table: run it after decontamination (--run_decontam true), " +
                  "after consensus (--start_from beginning), or use --start_from batch_correction with " +
                  "--decontam_otu_table and --metadata_file (audit C20)"
        }

        BatchCorrection(ch_batch_input, depsOf([params.batch_corr_script, "${params.scripts}/install_conqur.R"]))
    }

    // ============================================================================
    // COHORT QC SUMMARY
    // ============================================================================
    // One row per sample, assembled from whatever the requested stages produced. Runs
    // last so every contributing channel exists; a stage that did not run simply leaves
    // its columns NA.
    if (qcParts) {
        // Folded by hand, and seeded from the first element rather than
        // Channel.empty(), because of two parser constraints pulling opposite ways:
        // Nextflow 25+ rejects the spread operator `mix(*qcParts)` and `for` loops,
        // while 24.10 rejects `def x = Channel.empty()` inside the entry workflow
        // ("Variable `Channel` already defined in the process scope"). This form
        // compiles on both.
        def qcChannel = qcParts.first()
        qcParts.drop(1).each { part -> qcChannel = qcChannel.mix(part) }
        qcSummary(qcChannel.collect(), depsOf(["${params.scripts}/build_qc_summary.py"]))

        // The report reads the summary table, not the results tree, so the two cannot
        // disagree about what the run produced.
        if (!skipCohortReport) {
            // Optional inputs as value channels: [] when that stage did not run.
            if (!runValidation) { Channel.value([]).set { VALIDATION_MASK_REPORT } }
            if (!(params.start_from == 'beginning' && !skipClassification && !skipConsensus)) {
                Channel.value([]).set { CONSENSUS_STATUS }
            }
            cohortReport(qcSummary.out.summary, VALIDATION_MASK_REPORT, CONSENSUS_STATUS,
                         depsOf(["${params.scripts}/build_cohort_report.py"]))
        }
    }

    // ------------------------------------------------------------------------
    // Workflow completion handlers.
    // These live inside the entry workflow: the strict parser in Nextflow 25+
    // rejects statements at script level. Verified on 24.10.0 and 26.04.6 (manifest: >=24.10.0).
    // ------------------------------------------------------------------------
    // Inside the entry workflow, `workflow`, `params` and `projectDir` are all shadowed
    // within a nested closure and resolve to null. Bind them here and read the bindings
    // in the handlers.
    def wfMeta      = workflow
    def tracedirOut = params.tracedir.toString()
    def projectRoot = workflow.projectDir.toString()

    workflow.onComplete {
        log.info """
        ================================================================================
        Pipeline execution summary
        ================================================================================
        Completed at : ${wfMeta.complete}
        Duration     : ${wfMeta.duration}
        Success      : ${wfMeta.success}
        Exit status  : ${wfMeta.exitStatus}
        Error report : ${wfMeta.errorReport ?: '-'}
        ================================================================================
        """.stripIndent()

        // Provenance, completion half. Shelled out so the rendering and trace parsing
        // sit somewhere the test suite can reach. A provenance failure must never mask
        // the run's own outcome, hence the catch.
        try {
            def finalise = [
                'python3', "${projectRoot}/scripts/finalize_provenance.py",
                '--tracedir',    tracedirOut,
                '--session',     "${wfMeta.sessionId}_${wfMeta.runName}",
                '--success',     "${wfMeta.success}",
                '--exit-status', "${wfMeta.exitStatus}",
                '--duration',    "${wfMeta.duration}",
                '--error-message', (wfMeta.errorMessage ?: ''),
                // Authoritative counts straight from the session: onComplete can reach
                // trace.txt before the trace observer has flushed it, so the totals must
                // not depend on that file.
                '--succeeded', "${wfMeta.stats.succeededCount}",
                '--cached',    "${wfMeta.stats.cachedCount}",
                '--failed',    "${wfMeta.stats.failedCount}",
                '--ignored',   "${wfMeta.stats.ignoredCount}",
                // No --conda-cache-dir: the finaliser reads the envs this launch recorded
                // instead of every env in the cache (audit C06/D08).
            ]
            def proc = finalise.execute()
            proc.waitForProcessOutput(System.out, System.err)
        }
        catch (Exception problem) {
            System.err.println("WARN: could not finalise the provenance record: ${problem.message}")
        }
    }
    workflow.onError {
        log.info """
        ================================================================================
        Pipeline execution error
        ================================================================================
        ${wfMeta.errorMessage}
        ================================================================================
        """.stripIndent()
    }
}
