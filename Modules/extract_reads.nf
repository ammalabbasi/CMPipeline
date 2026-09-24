nextflow.enable.dsl=2

process extractReads {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    tag "$sampleID"
    scratch true
    label 'extract_reads'
    // Intermediate read/alignment files are published only with --save_intermediates true (audit C11);
    publishDir("${params.unmapped_bam_dir}", mode: 'copy', saveAs: { fn -> (params.save_intermediates.toString() == 'true' || !(fn =~ /\.(fastq|fq|fasta|fa|sam|bam|bowtie2\.bz2)(\.gz|\.bz2)?$/)) ? fn : null })
    conda "${params.samtools_env}"

    input:
    tuple val(sampleID), val(inputType), path(alignment), path(reference_files)

    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple val(sampleID), path("${sampleID}.UNMAPPED.fastq.gz"), emit: reads
    // Library denominator for BAM/CRAM (audit C07 / pksProfiler F09): primary records in the input
    // and how many of them were unmapped and extracted. Counts only, no read content.
    path("${sampleID}.extract_counts.tsv"), emit: counts

    script:
    def reference = reference_files ? reference_files[0] : null
    """
    set -euo pipefail

    READS="${sampleID}.UNMAPPED.fastq.gz"
    FORMAT="\$(htsfile "${alignment}" 2>/dev/null || true)"
    case "\$FORMAT" in
        *CRAM*) DETECTED=cram ;;
        *BAM*)  DETECTED=bam ;;
        *) echo "ERROR: ${sampleID}: input is neither BAM nor CRAM according to htsfile: \$FORMAT" >&2; exit 1 ;;
    esac

    if [[ "${inputType}" == "bam" && "\$DETECTED" != "bam" ]]; then
        echo "ERROR: ${sampleID}: BAM was requested but the input is CRAM" >&2
        exit 1
    fi
    if [[ "${inputType}" == "cram" && "\$DETECTED" != "cram" ]]; then
        echo "ERROR: ${sampleID}: CRAM was requested but the input is BAM" >&2
        exit 1
    fi

    REFERENCE_ARGS=()
    if [[ "\$DETECTED" == "cram" && -n "${reference ?: ''}" ]]; then
        python "${params.scripts}/validate_cram_reference.py" \
            --alignment "${alignment}" --reference "${reference ?: ''}"
        REFERENCE_ARGS=(-T "${reference ?: ''}")
    fi

    samtools quickcheck -v "${alignment}"
    # `samtools fastq -N -o FILE` writes only READ1/READ2 to FILE and sends category-0
    # records -- reads flagged neither READ1 nor READ2 -- to STDOUT, which in a Nextflow
    # task is .command.out. That silently dropped those reads from the FASTQ *and* put
    # read sequences in the task log. Verified on samtools 1.21 with a synthetic BAM and
    # CRAM: `-o FILE` wrote 3 of 5 reads and leaked the other 2; naming neither -o nor -0
    # writes all 5 and leaks nothing. Same defect as pksProfiler F01.
    # The locked env (conda_envs/lock/samtools_env) is samtools 1.22.1 / htslib 1.24; the routing is
    # re-checked on exactly that build by cmpipeline/CHECK_ROUND2.sh (audit D11), and
    # tests/test_repo_hygiene.py fails if the lockfile's samtools moves without this line.
    # Stream every category to stdout and compress it, so nothing is routed away.
    # Primary records are counted in the SAME decode pass (tee into two counting fifos), so the
    # library denominator costs no second read of a 50-90 GB CRAM (audit C07).
    mkfifo primary.fifo unmapped.fifo
    samtools view -c primary.fifo  > primary.count  & P_PRIMARY=\$!
    samtools view -c unmapped.fifo > unmapped.count & P_UNMAPPED=\$!
    if ! samtools view -@ "${task.cpus}" -F 2304 -u \
        "\${REFERENCE_ARGS[@]}" "${alignment}" |
        tee primary.fifo |
        samtools view -f 4 -u - |
        tee unmapped.fifo |
        samtools fastq -@ "${task.cpus}" -N - |
        bgzip -@ "${task.cpus}" -c > "\$READS"; then
        echo "ERROR: Could not decode ${alignment}. For CRAM, provide the matching --cram_reference if required." >&2
        exit 1
    fi
    wait \$P_PRIMARY
    wait \$P_UNMAPPED
    printf 'sample\tlibrary_primary_records\textracted_unmapped_records\n%s\t%s\t%s\n' \
        "${sampleID}" "\$(< primary.count)" "\$(< unmapped.count)" > "${sampleID}.extract_counts.tsv"
    rm -f primary.fifo unmapped.fifo primary.count unmapped.count
    gzip -t "\$READS"
    """
}
