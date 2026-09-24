nextflow.enable.dsl=2

/*
 * Cohort QC report.
 *
 * The previous version of this module could not have worked: it passed --no-report to
 * MultiQC while declaring `multiqc_report.html` as an output, and declared
 * `multiqc_data.tsv` where --data-dir writes a `multiqc_data/` directory. It was also
 * never imported by main.nf, so neither fault ever surfaced.
 *
 * Inputs are staged into numbered subdirectories: FASTQC runs five times per sample
 * (raw, filtered, hg38, t2t, pangenome) and its zip is named after the input file, so a
 * flat staging can collide. MultiQC scans recursively.
 */
process multiqc {

    label 'process_medium'
    publishDir("${params.multiqc_dir}", mode: 'copy')
    conda "${params.multiqc_env}"

    input:
    path(qc_files, stageAs: "qc_?/*")

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"),        emit: data

    script:
    """
    set -euo pipefail
    multiqc -f -p \\
        --data-dir --data-format tsv \\
        --cl-config "max_table_rows: 10000" \\
        --outdir ./ .
    """
}
