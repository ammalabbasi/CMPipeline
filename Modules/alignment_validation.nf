// =============================================================================
// ALIGNMENT VALIDATION (opt-in: --run_alignment_validation)
// =============================================================================
//
// Consensus asks whether two classifiers agreed; decontam asks whether the reagents put a taxon
// there. This asks whether a taxon's reads actually come from its genome, and catches
// misassignment: rRNA/plasmid pile-ups, reads of an absent relative, and "shadow" species that
// only share reads with a real one. Design: DESIGN_alignment_validation.md (v2, 2026-09-23);
// thresholds from the simulation benchmark in alignment_benchmark/RESULTS.md.
//
//   validationSpeciesList -> buildValidationRef (cached per species list)
//     -> alignForValidation (batches of samples, one minimap2 run each) -> applyValidationMask
//
// Bracken counts stay the quantity downstream; validation only zeroes (sample, species) cells
// that FAILED. minimap2 counts are reported, never used.
// =============================================================================

nextflow.enable.dsl=2

// RefSeq snapshot month for the validation reference: --validation_refseq_snapshot, else the month
// this run started (decision 26). Both processes call it, so they always agree within a run.
def validationSnapshot() {
    params.validation_refseq_snapshot ?: workflow.start.toString().substring(0, 7)
}

process validationSpeciesList {
    label 'process_low'
    conda "${params.validation_env}"
    publishDir "${params.validation_dir}", mode: 'copy'

    input:
    path species_table
    path taxdb

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple env('REF_KEY'), path("validation_species.tsv"), emit: species

    script:
    // Everything that determines the reference goes into its storeDir cache key: the species
    // (taxids, in the script), the build settings, and the builder/helper scripts' digests (deps).
    // A changed builder or setting must not silently reuse an old reference (audit B04/C03).
    // The RefSeq snapshot MONTH is in the key too (decision 26): within a month the cache is
    // reused; a new month gets a new key and a fresh panel instead of silently replacing the old
    // one. Per-genome summary dates are in validation_reference_manifest.tsv.
    def settings = "max_genomes=${params.validation_max_genomes};ani=${params.validation_dedup_ani};" +
                   "max_missing=${params.validation_max_missing_frac};t2t=GCF_009914755.1;" +
                   "refseq=${validationSnapshot()};deps=${deps};v4"
    """
    python "${params.scripts}/validation_species_list.py" \\
        --species_table "${species_table}" \\
        --taxdb "${taxdb}" \\
        --min_support ${params.validation_min_support} \\
        --key_settings "${settings}"
    REF_KEY=\$(cat ref_key.txt)
    """
}

process buildValidationRef {
    tag "ref-${ref_key}"
    cpus 16
    memory { 96.GB * task.attempt }
    time { 12.h * task.attempt }
    conda "${params.validation_env}"
    // Cached per species list + build settings: a re-run, or another run with the same
    // species, reuses the index without rebuilding or re-downloading.
    storeDir { "${params.validation_cache_dir}/ref-${ref_key}" }

    input:
    tuple val(ref_key), path(species)

    output:
    tuple val(ref_key),
          path("validation_reference.mmi"),
          path("contig_labels.tsv"),
          path("validation_reference_manifest.tsv"),
          path("missing_species.tsv"), emit: ref

    script:
    """
    set -euo pipefail
    python "${params.scripts}/build_validation_ref.py" \\
        --species "${species}" \\
        --cache_dir "${params.validation_cache_dir}/downloads" \\
        --max_genomes ${params.validation_max_genomes} \\
        --dedup_ani ${params.validation_dedup_ani} \\
        --max_missing_frac ${params.validation_max_missing_frac} \\
        --snapshot ${validationSnapshot()} \\
        --max_ref_gb ${params.validation_max_ref_gb} \\
        --threads ${task.cpus}

    # ONE-PART index. minimap2's default batch (-I 8G) split the benchmark reference in two,
    # mapped every read against each part separately and reported it twice: microbes and human
    # never competed. -I is set above the reference size, and a split is refused.
    GB=\$(( \$(stat -c %s validation_reference.fna) / 1000000000 + 2 ))
    minimap2 -x sr -I \${GB}G -t ${task.cpus} -d validation_reference.mmi validation_reference.fna 2> index.log
    PARTS=\$(grep -c "loaded/built the index" index.log || true)
    if [ "\$PARTS" -ne 1 ]; then
        echo "ERROR: validation index has \$PARTS parts (expected 1); competition would be broken" >&2
        exit 1
    fi
    rm -f validation_reference.fna
    """
}

process alignForValidation {
    tag "${sample_ids.size()} sample(s)"
    cpus 16
    // Scales with the reference (audit B05). Benchmarked peak RSS = index size + ~9 GB (12 GB index ->
    // 21 GB; 24 GB -> 32 GB), so request 1.2 x index + 12 GB, doubled per retry.
    // No Groovy (long) cast: the Nextflow 26 strict parser rejects it (caught by tests/integration).
    memory { Math.ceil(index.size() / 1e9 * 1.2 + 12).toLong().GB * task.attempt }
    time { 12.h * task.attempt }
    conda "${params.validation_env}"
    scratch true
    publishDir "${params.validation_dir}/per_sample", mode: 'copy'

    input:
    tuple val(sample_ids), path(reads)
    tuple val(ref_key), path(index), path(contig_labels), path(manifest), path(missing)
    path species_table
    path groups

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    path("*.validation_stats.tsv"), emit: stats
    path("*.validated_reads.fasta.gz"), optional: true, emit: reads

    script:
    // One minimap2 run per batch: the index (tens of GB) loads once. Read names are prefixed
    // <sample>:: so the scorer can split the stream back into samples.
    def pairs = [sample_ids, reads].transpose()
    def stream = pairs.collect { id, f ->
        "gzip -dc \"${f}\" | awk -v S='${id}' 'NR%4==1{sub(/^@/, \"@\" S \"::\")} {print}'"
    }.join('\n        ')
    def groups_arg = groups ? "--groups \"${groups}\"" : ""
    def export_arg = params.validation_export_reads.toString() == 'true' ? "--export_reads" : ""
    """
    set -euo pipefail
    {
        ${stream}
    } | minimap2 -ax sr -t ${task.cpus} -p 0.99 -N 100 --secondary=yes "${index}" /dev/stdin 2> minimap2.log \\
      | python "${params.scripts}/validation_score.py" \\
            --samples "${sample_ids.join(',')}" \\
            --contigs "${contig_labels}" \\
            --species_table "${species_table}" \\
            ${groups_arg} \\
            --min_support ${params.validation_min_support} \\
            --min_confirmed ${params.validation_min_confirmed} \\
            --min_bins_ratio ${params.validation_min_bins_ratio} \\
            --min_identity ${params.validation_min_identity} \\
            --min_unique_frac ${params.validation_min_unique_frac} \\
            ${export_arg}
    """
}

process applyValidationMask {
    label 'process_low'
    conda "${params.validation_env}"
    publishDir "${params.validation_dir}", mode: 'copy'

    input:
    path species_table
    path stats
    path groups

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    path("validated.species.tsv"), emit: species
    path("validated.genus.tsv"), emit: genus
    path("validation_matrix.tsv")
    path("validation_stats.all.tsv")
    path("validation_mask_report.tsv"), emit: report

    script:
    def groups_arg = groups ? "--groups \"${groups}\"" : ""
    """
    python "${params.scripts}/apply_validation_mask.py" \\
        --species_table "${species_table}" \\
        --stats ${stats} \\
        ${groups_arg}
    """
}

process compareDecontamLevels {
    label 'process_low'
    conda "${params.validation_env}"
    publishDir "${params.validation_dir}", mode: 'copy'

    input:
    path species_final
    path genus_final

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    path("decontam_species_vs_genus.tsv")

    script:
    """
    python "${params.scripts}/compare_decontam_levels.py" \\
        --species_final "${species_final}" --genus_final "${genus_final}"
    """
}
