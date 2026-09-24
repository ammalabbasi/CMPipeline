nextflow.enable.dsl=2

process mapReads {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    tag "$sampleID"
    scratch true
    label 'mapReads'
    // Intermediate read/alignment files are published only with --save_intermediates true (audit C11);
    // the final host-depleted reads are always published.
    publishDir("${params.mapped_reads_dir}", mode: 'copy', saveAs: { fn -> (params.save_intermediates.toString() == 'true' || fn.endsWith('.host_depleted.fastq.gz') || !(fn =~ /\.(fastq|fq|fasta|fa|sam|bam|bowtie2\.bz2)(\.gz|\.bz2)?$/)) ? fn : null })
    conda "${params.minimap2_env}"

    input:
    tuple val(sampleID), path(reads_fastq)

    // Named emits: a second output would turn the process handle into a multi-channel
    // object and break `mapReads(...).set{}` at the call site.
    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple val(sampleID),
        path("${sampleID}.hg38.fastq.gz"),
        path("${sampleID}.hg38.t2t.fastq.gz"),
        path("${sampleID}.host_depleted.fastq.gz"), emit: reads
    // Read retention at each depletion stage. This is the pipeline's central QC claim
    // and until now it was not recorded anywhere.
    path("${sampleID}.depletion_counts.tsv"), emit: counts

    script:
    """
    set -euo pipefail
    gzip -t "${reads_fastq}"

    minimap2 -2 -ax sr -t "${task.cpus}" "${params.hg38_db}" "${reads_fastq}" |
        samtools fastq -@ "${task.cpus}" -f 4 -F 2304 | gzip -c > "${sampleID}.hg38.fastq.gz"

    minimap2 -2 -ax sr -t "${task.cpus}" "${params.t2t_phix_db}" "${sampleID}.hg38.fastq.gz" |
        samtools fastq -@ "${task.cpus}" -f 4 -F 2304 | gzip -c > "${sampleID}.hg38.t2t.fastq.gz"

    if [[ -n "${params.pangenome_db ?: ''}" ]]; then
        if [[ -d "${params.pangenome_db}" ]]; then
            mapfile -t INDEXES < <(find "${params.pangenome_db}" -type f -name '*.mmi' -print | sort)
        elif [[ -f "${params.pangenome_db}" ]]; then
            INDEXES=("${params.pangenome_db}")
        else
            echo "ERROR: Pangenome index path does not exist: ${params.pangenome_db}" >&2; exit 1
        fi
        [[ "\${#INDEXES[@]}" -gt 0 ]] || { echo "ERROR: No .mmi files found under ${params.pangenome_db}" >&2; exit 1; }
        cp "${sampleID}.hg38.t2t.fastq.gz" "${sampleID}.current.fastq.gz"
        for mmi in "\${INDEXES[@]}"; do
            minimap2 -2 -ax sr -t "${task.cpus}" "\$mmi" "${sampleID}.current.fastq.gz" |
                samtools fastq -@ "${task.cpus}" -f 4 -F 2304 | gzip -c > "${sampleID}.next.fastq.gz"
            mv "${sampleID}.next.fastq.gz" "${sampleID}.current.fastq.gz"
        done
        mv "${sampleID}.current.fastq.gz" "${sampleID}.host_depleted.fastq.gz"
    else
        cp "${sampleID}.hg38.t2t.fastq.gz" "${sampleID}.host_depleted.fastq.gz"
    fi

    # Read retention per stage: four integers per sample, no sequence written anywhere.
    # It is a second decompression pass over files that were just written, so it is
    # cheap next to minimap2 and is the same pattern the module already uses for
    # `gzip -t`. Its real cost will be visible per process in pipeline_info/trace.txt.
    count_reads() {
        gzip -dc "\$1" | awk 'END { print int(NR / 4) }'
    }
    {
        printf 'sample\\tstage\\treads\\n'
        printf '%s\\tinput\\t%s\\n'          "${sampleID}" "\$(count_reads "${reads_fastq}")"
        printf '%s\\tafter_hg38\\t%s\\n'     "${sampleID}" "\$(count_reads "${sampleID}.hg38.fastq.gz")"
        printf '%s\\tafter_t2t_phix\\t%s\\n' "${sampleID}" "\$(count_reads "${sampleID}.hg38.t2t.fastq.gz")"
        printf '%s\\thost_depleted\\t%s\\n'  "${sampleID}" "\$(count_reads "${sampleID}.host_depleted.fastq.gz")"
    } > "${sampleID}.depletion_counts.tsv"
    """
}
