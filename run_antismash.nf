#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Standalone antiSMASH workflow for metagenome assemblies and/or MAGs
 *
 * Two modes:
 *   contigs - Run on MEGAHIT assembly contigs directly (recommended for
 *             cancer microbiome data where binning often fails)
 *   mags    - Run on binned MAGs from nf-core/mag GenomeBinning output
 *
 * Usage:
 *   nextflow run run_antismash.nf --mode contigs -profile conda -resume
 *   nextflow run run_antismash.nf --mode mags    -profile conda -resume
 */

// ── Parameters ─────────────────────────────────────────────────────────────

// Mode: 'contigs' or 'mags'
params.mode            = "contigs"

// Path to nf-core/mag output
params.mag_output_dir  = "${projectDir}/RESULTS/MAG"

// Output
params.antismash_dir   = "${projectDir}/RESULTS/ANTISMASH"

// Conda env
params.antismash_env   = "${projectDir}/conda_envs/antismash_env.yml"

// ── Include module ─────────────────────────────────────────────────────────
include { antiSMASH } from './Modules/antismash.nf'

// ── Workflow ───────────────────────────────────────────────────────────────
workflow {

    log.info """
    ===================================================
    antiSMASH — mode: ${params.mode}
    ===================================================
    MAG output dir : ${params.mag_output_dir}
    Output dir     : ${params.antismash_dir}
    ===================================================
    """.stripIndent()

    if (params.mode == "contigs") {
        /*
         * Run on MEGAHIT assembly contigs directly.
         * Files: RESULTS/MAG/Assembly/MEGAHIT/MEGAHIT-{sample}.contigs.fa.gz
         * ID extracted as sample name (e.g. PD44724a from MEGAHIT-PD44724a).
         */
        input_fastas = Channel.fromPath(
            "${params.mag_output_dir}/Assembly/MEGAHIT/MEGAHIT-*.contigs.fa.gz",
            checkIfExists: true
        )
        .map { fasta ->
            // MEGAHIT-PD44724a.contigs.fa.gz -> PD44724a
            def sample_id = fasta.name.replaceAll(/^MEGAHIT-/, '').replaceAll(/\.contigs\.fa\.gz$/, '')
            tuple(sample_id, fasta)
        }
    } else if (params.mode == "mags") {
        /*
         * Run on binned MAGs from nf-core/mag GenomeBinning output.
         * Searches DASTool, MetaBAT2, CONCOCT bin directories.
         */
        input_fastas = Channel.fromPath(
            [
                "${params.mag_output_dir}/**/GenomeBinning/**/bins/*.fa.gz",
                "${params.mag_output_dir}/**/GenomeBinning/**/bins/*.fa",
                "${params.mag_output_dir}/**/GenomeBinning/**/bins/*.fasta",
                "${params.mag_output_dir}/**/GenomeBinning/**/bins/*.fasta.gz",
            ],
            checkIfExists: false
        )
        .map { fasta ->
            def mag_id = fasta.baseName.replaceAll(/\.fa(sta)?$/, '').replaceAll(/\.gz$/, '')
            tuple(mag_id, fasta)
        }
    } else {
        error "Unknown mode '${params.mode}'. Use --mode contigs or --mode mags"
    }

    input_fastas.view { id, fasta -> "Found ${params.mode == 'contigs' ? 'assembly' : 'MAG'}: ${id} -> ${fasta}" }

    // Run antiSMASH
    antiSMASH(input_fastas)

    // Summarize
    antiSMASH.out.json
        .map { id, json_files -> id }
        .collect()
        .view { ids -> "\nProcessed ${ids.size()} ${params.mode == 'contigs' ? 'assemblies' : 'MAGs'} through antiSMASH" }
}
