// =============================================================================
// BATCH CORRECTION MODULE
// =============================================================================
//
// Batch effect correction and normalization of microbiome data using ConQuR
//
// =============================================================================

nextflow.enable.dsl=2

process BatchCorrection {

    tag "${prefix}"
    label 'process_high'

    conda "${params.batch_corr_env}"
    // Dynamic path: see the note in decontamination.nf.
    publishDir path: { "${params.batch_corr_dir}/${prefix}" }, mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata)

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple val(prefix),
          path("corrected/*"),
          path("normalized/*.tsv"),
          path("pcoa_plots/*"),
          path("permanova_summary.tsv"),
          path("batch_correction_sample_status.tsv"), emit: corrected

    script:
    // `--tumor_only false` arrives as the String "false" on Nextflow 25+, which is
    // truthy in Groovy. Compare explicitly rather than testing truthiness.
    def tumor_only_flag = (params.tumor_only.toString().trim().toLowerCase() in ['true', 'yes', 'on', '1']) ? "--tumor_only" : ""
    def method_flag = params.batch_corr_method ?: "tune"
    def allow_flag = params.batch_corr_allow_uncorrected.toString() == 'true' ? "--allow_uncorrected" : ""
    """
    # ConQuR is installed from GitHub, not from the conda YAML, so an environment built
    # only with `conda env create` is missing it. Fail here with the fix rather than
    # partway through the R script.
    # Presence is not enough: it must be the pinned fork the script was written against (audit A14).
    Rscript -e 'if (!requireNamespace("ConQuR", quietly = TRUE)) stop("ConQuR is not installed in this environment (env: ", Sys.getenv("CONDA_PREFIX"), "). Run: Rscript ${params.scripts}/install_conqur.R with the batch-correction env active"); sha <- utils::packageDescription("ConQuR")\$RemoteSha; if (is.null(sha) || !startsWith(sha, "${params.conqur_sha}")) stop("ConQuR in this environment is not ivartb/ConQuR_par@${params.conqur_sha} (found: ", if (is.null(sha)) "no RemoteSha" else sha, "; env: ", Sys.getenv("CONDA_PREFIX"), "). Run: CONQUR_FORCE=1 Rscript ${params.scripts}/install_conqur.R")'

    # Every value quoted: the TCGA preset's type column is 'Sample Type_x' (audit A03).
    Rscript "${params.batch_corr_script}" \\
        --otu "${otu_table}" \\
        --meta "${metadata}" \\
        --prefix "${prefix}" \\
        --batch_column "${params.batch_column}" \\
        --covariates "${params.batch_corr_covariates}" \\
        --type_column "${params.type_column}" \\
        --sample_id_column "${params.sample_id_column}" \\
        --taxon_column "${params.taxon_column}" \\
        --tumor_value "${params.batch_corr_tumor_values}" \\
        ${allow_flag} \\
        ${tumor_only_flag} \\
        --phase ${params.phase} \\
        --r2_threshold ${params.r2_threshold} \\
        --method ${method_flag}
    """
}
