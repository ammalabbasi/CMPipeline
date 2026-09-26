// =============================================================================
// BATCH CORRECTION MODULE
// =============================================================================
//
// Batch effect correction and normalization of microbiome data using ConQuR
//
// =============================================================================

nextflow.enable.dsl=2

process BatchCorrection {

    tag "${prefix}_${rank}"
    label 'process_high'

    conda "${params.batch_corr_env}"
    // Dynamic path: see the note in decontamination.nf. One directory per rank, e.g.
    // decontam_run_species/ and decontam_run_genus/. The QC-only ledger copy is not published.
    publishDir path: { "${params.batch_corr_dir}/${prefix}_${rank}" }, mode: 'copy', overwrite: true,
               saveAs: { name -> name.endsWith('.batch_sample_status.tsv') ? null : name }

    input:
    // rank: species | genus -- the level of THIS table.
    tuple val(prefix), path(otu_table), path(metadata), val(rank)

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
          path("batch_correction_sample_status.tsv"),
          path("batch_correction_taxon_status.tsv"),
          path("batch_correction_quality.tsv"),
          path("batch_correction_methods.tsv"),
          path("batch_correction_run_settings.tsv"), emit: corrected
    path("${rank}.batch_sample_status.tsv"), emit: qc_status

    script:
    // `--tumor_only false` arrives as the String "false" on Nextflow 25+, which is
    // truthy in Groovy. Compare explicitly rather than testing truthiness.
    def tumor_only_flag = (params.tumor_only.toString().trim().toLowerCase() in ['true', 'yes', 'on', '1']) ? "--tumor_only" : ""
    def method_flag = params.batch_corr_method ?: "tune"
    def allow_flag = params.batch_corr_allow_uncorrected.toString() == 'true' ? "--allow_uncorrected" : ""
    """
    export R_LIBS="\$CONDA_PREFIX/lib/R/library"
    export R_LIBS_USER="\$CONDA_PREFIX/lib/R/library"
    export R_LIBS_SITE="\$CONDA_PREFIX/lib/R/library"
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
    # ConQuR comes from the repository's local conda channel (conda_envs/local_channel/). A
    # pre-built env passed with --batch_corr_env may lack it: fail here with the fix rather
    # than partway through the R script.
    # Presence is not enough: it must be the pinned fork the script was written against (audit A14).
    Rscript --vanilla -e 'if (!requireNamespace("ConQuR", quietly = TRUE)) stop("ConQuR is not installed in this environment (env: ", Sys.getenv("CONDA_PREFIX"), "). Use the default environment (conda_envs/lock/batch_correction_env.linux-64.txt, built by Nextflow) or rebuild the prefix with scripts/setup_batch_env.sh"); sha <- utils::packageDescription("ConQuR")\$RemoteSha; if (is.null(sha) || !startsWith(sha, "${params.conqur_sha}")) stop("ConQuR in this environment is not ivartb/ConQuR_par@${params.conqur_sha} (found: ", if (is.null(sha)) "no RemoteSha" else sha, "; env: ", Sys.getenv("CONDA_PREFIX"), "). Use the default environment or rebuild the prefix with scripts/setup_batch_env.sh")'

    # Every value quoted: the TCGA preset's type column is 'Sample Type_x' (audit A03).
    Rscript --vanilla "${params.batch_corr_script}" \\
        --otu "${otu_table}" \\
        --meta "${metadata}" \\
        --prefix "${prefix}_${rank}" \\
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
        --method ${method_flag} \\
        --seed ${params.batch_corr_seed} \\
        --threads ${task.cpus}

    # Rank-named copy of the sample ledger for the cohort QC summary.
    cp batch_correction_sample_status.tsv "${rank}.batch_sample_status.tsv"
    """
}
