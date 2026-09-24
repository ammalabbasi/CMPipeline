// =============================================================================
// DECONTAMINATION MODULE
// =============================================================================
//
// Statistical identification and removal of contaminants from microbiome data
// using the decontam R package
//
// =============================================================================

nextflow.enable.dsl=2

process Decontamination {

    tag "${prefix}"
    label 'process_medium'

    conda "${params.decontam_env}"
    // Dynamic path: `prefix` is an input value, so the path must be a closure.
    // Nextflow 25+ evaluates a plain GString directive eagerly and fails with
    // "No such variable: prefix".
    publishDir path: { "${params.decontam_dir}/${prefix}" }, mode: 'copy', overwrite: true

    input:
    // library_sizes: per-sample Bracken totals from consensus_taxa, or [] when unavailable
    // (--skip_consensus, or --start_from decontam without --decontam_library_sizes).
    // rank: genus | species -- the level of THIS table. With alignment validation on, decontam
    // runs once per level on the validated tables; otherwise it is --taxonomic_rank.
    // library_sizes may be the very same file as otu_table (consensus skipped), so it is staged
    // under its own folder to avoid a name collision.
    tuple val(prefix), path(otu_table), path(metadata), path(library_sizes, stageAs: 'library/*'), val(rank)

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple val(prefix),
          path("${prefix}.decontam_pkg_decontaminated.csv"),
          path("${prefix}.validated_metadata.tsv"),
          path("${prefix}.sample_validation.tsv"),
          path("${prefix}.batch_eligibility.tsv"),
          path("${prefix}.contaminants_by_batch.tsv"),
          path("${prefix}.contaminants_final.tsv"),
          path("${prefix}.filter_report.tsv"),
          path("${prefix}.run_manifest.tsv"),
          path("decontamination_plots/*"),
          path("${prefix}.decontamination_summary.txt"), emit: decontaminated
    tuple val(prefix),
          path("${prefix}.decontam_pkg_decontaminated.csv"),
          path("${prefix}.validated_metadata.tsv"), emit: for_batch_correction

    script:
    def library_arg = library_sizes ? "--library_table \"${library_sizes}\"" : ""
    def unknown_arg = params.decontam_allow_unknown_library.toString() == 'true' ? "--allow_unknown_library" : ""
    """
    set -euo pipefail
    Rscript "${params.decontam_script}" \\
        --otu_table "${otu_table}" \\
        --metadata "${metadata}" \\
        --prefix "${prefix}" \\
        --start_mode "${params.start_from}" \\
        --environment_spec "${params.decontam_env}" \\
        --nextflow_version "${workflow.nextflow.version}" \\
        --taxon_column "${params.taxon_column}" \\
        --taxonomic_rank "${rank}" \\
        --sample_id_column "${params.sample_id_column}" \\
        --batch_column "${params.batch_column}" \\
        --type_column "${params.type_column}" \\
        --control_role "${params.control_role}" \\
        --positive_values "${params.positive_values}" \\
        --control_values "${params.control_values}" \\
        --unselected_type_policy "${params.unselected_type_policy}" \\
        --threshold "${params.decontam_threshold}" \\
        --min_prevalence "${params.decontam_min_prevalence}" \\
        --min_abundance "${params.decontam_min_abundance}" \\
        --min_batches "${params.decontam_min_batches}" \\
        --min_total_per_batch "${params.min_total_per_batch}" \\
        --min_positive_per_batch "${params.min_positive_per_batch}" \\
        --min_control_per_batch "${params.min_control_per_batch}" \\
        --zero_total_policy "${params.decontam_zero_total_policy}" \\
        --unevaluable_policy "${params.decontam_unevaluable_policy}" \\
        --min_library ${params.decontam_min_library} ${unknown_arg} ${library_arg}
    """
}
