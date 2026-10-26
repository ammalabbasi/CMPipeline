nextflow.enable.dsl=2

process extractReads {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    tag "$sampleID"
    scratch true
    label 'extract_reads'
    // Intermediate read/alignment files are published only with --save_intermediates true (audit C11);
    // the rescued-read name list (--extract_low_score_mapped) is an intermediate too: read names only.
    publishDir("${params.unmapped_bam_dir}", mode: 'copy', saveAs: { fn -> (params.save_intermediates.toString() == 'true' || !(fn =~ /\.(fastq|fq|fasta|fa|sam|bam|bowtie2\.bz2)(\.gz|\.bz2)?$/ || fn.endsWith('.low_as_read_ids.txt.gz'))) ? fn : null })
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
    // and how many of them were unmapped and extracted, plus the low-AS reads the rescue added.
    // Counts only, no read content.
    path("${sampleID}.extract_counts.tsv"), emit: counts
    // --extract_low_score_mapped: the upstream aligner, its version and bwa -T from @PG (never the
    // CL itself), and the mode this sample got and why. One row, metadata only.
    path("${sampleID}.upstream_aligner.tsv"), emit: aligner
    // Input provenance (audit X13, Ludmil E3): @SQ fingerprint and contig classes (decoy / viral /
    // HLA / ALT, chrEBV), the @PG chain as program:version, bwa mem options reduced to letters and
    // numbers, and any upstream step that removed unmapped reads. From the same header text as the
    // rescue decision; never the CL, a read group, a path or a sample barcode.
    path("${sampleID}.input_provenance.tsv"), emit: provenance
    // Names of the rescued low-AS reads exactly as `samtools fastq -N` wrote them, so lowAsSurvival
    // can count them after host depletion and classification. Only when the rescue is on.
    tuple val(sampleID), path("${sampleID}.low_as_read_ids.txt.gz"), emit: low_as_ids, optional: true
    // QC floor (pksProfiler COHORT_CHECKS A3): pass, or why this input cannot be assessed --
    // unmapped_filtered_upstream (the @PG chain shows a step that removed unmapped reads) or
    // no_unmapped_records (zero FLAG 4 records). main.nf stops such a sample here
    // (--qc_floor_action exclude), so an empty profile is never reported as a microbial negative.
    tuple val(sampleID), env('INPUT_QC_FLOOR'), emit: floor

    script:
    def reference = reference_files ? reference_files[0] : null
    """
    set -euo pipefail
    ${ToolVersions.record(task, params, 'samtools htsfile:htslib bgzip:htslib python')}

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

    # Low-score mapped rescue (--extract_low_score_mapped auto|on|off; Ludmil 1 Oct, L8). bwa mem
    # reports a read unmapped only when its best local score is below -T (default 30). GDC aligned
    # TCGA with `bwa mem -T 0`, so a microbial read with a chance ~19 bp human match is reported
    # MAPPED with AS < 30 and FLAG 4 extraction never sees it. `auto` turns the rescue on only when
    # the header's @PG runs bwa mem / bwa-mem2 mem with -T < 30; STAR, minimap2, bwa aln, bwa mem at
    # the default -T, or no @PG leave it off. Decided per sample from the header alone (the CL is
    # read in memory and never written; scripts/low_score_rescue.py). Ported from pksProfiler E1.
    # The header is read from the alignment ONCE: the same text also gives the input provenance row
    # (--provenance; scripts/header_provenance.py, audit X13). --no-PG: samtools would otherwise add
    # its own @PG for this header dump, and the pipeline's samtools would appear in the input's chain.
    samtools view -H --no-PG "${alignment}" > input.header.sam
    LOW_AS="\$(python "${params.scripts}/low_score_rescue.py" decide \
        --header input.header.sam --sample "${sampleID}" \
        --mode "${params.extract_low_score_mapped}" \
        --output "${sampleID}.upstream_aligner.tsv" \
        --provenance "${sampleID}.input_provenance.tsv")"
    rm -f input.header.sam

    # `samtools fastq -N -o FILE` writes only READ1/READ2 to FILE and sends category-0
    # records -- reads flagged neither READ1 nor READ2 -- to STDOUT, which in a Nextflow
    # task is .command.out. That silently dropped those reads from the FASTQ *and* put
    # read sequences in the task log. Verified on samtools 1.21 with a synthetic BAM and
    # CRAM: `-o FILE` wrote 3 of 5 reads and leaked the other 2; naming neither -o nor -0
    # writes all 5 and leaks nothing. Same defect as pksProfiler F01.
    # The locked env (conda_envs/lock/samtools_env) is samtools 1.24 / htslib 1.24; the routing is
    # re-checked on exactly that build, for both the FLAG 4 and the AS < 30 branch, with BAM and CRAM, by
    # audits/2026-10-26_toolbump/D11_samtools_1.24_routing.sh (audit D11; 4/4 PASS, 2026-10-26), and
    # tests/test_repo_hygiene.py fails if the lockfile's samtools moves without this line.
    # Stream every category to stdout and compress it, so nothing is routed away.
    # Primary records are counted in the SAME decode pass (tee into counting fifos), so the
    # library denominator costs no second read of a 50-90 GB CRAM (audit C07).
    mkfifo primary.fifo
    samtools view -c primary.fifo  > primary.count  & P_PRIMARY=\$!
    LOW_AS_READS=0
    NO_AS_RECORDS=NA
    if [[ "\$LOW_AS" == "on" ]]; then
        # Same single decode. The second samtools view passes on FLAG 4 records, records with
        # AS:i < 30, and mapped records with no AS at all (only to be counted: the awk drops them).
        # -F 2304 on the decode already removed secondary and supplementary records, so only
        # PRIMARY mapped records can be rescued. The awk counts each class and writes each rescued
        # read's name with the /1 or /2 that samtools fastq -N gives it (READ1 xor READ2; neither
        # or both: no suffix). samtools fastq then writes names, mate suffixes and reverse-strand
        # orientation exactly as in the off branch. A fifo counts the FASTQ reads for an add-up check.
        mkfifo fastq.fifo
        awk 'END { print int(NR / 4) }' < fastq.fifo > fastq.count & P_FASTQ=\$!
        if ! samtools view -@ "${task.cpus}" -F 2304 -u \
            "\${REFERENCE_ARGS[@]}" "${alignment}" |
            tee primary.fifo |
            samtools view -h -e '(flag & 4) || !exists([AS]) || [AS] < 30' - |
            awk -F '\\t' -v min_as=30 -v counts=low_as.classes -v ids=low_as.ids '
                /^@/ { print; next }
                {
                    if (int(\$2 / 4) % 2 == 1) { unmapped++; print; next }
                    score = ""
                    for (i = 12; i <= NF; i++) if (\$i ~ /^AS:[if]:/) { score = substr(\$i, 6); break }
                    if (score == "") { no_as++; next }
                    if (score + 0 < min_as + 0) {
                        low++
                        r1 = int(\$2 / 64) % 2; r2 = int(\$2 / 128) % 2; sfx = ""
                        if (r1 && !r2) sfx = "/1"; else if (r2 && !r1) sfx = "/2"
                        print \$1 sfx > ids
                        print; next
                    }
                    other++
                }
                END { printf "%d\\t%d\\t%d\\t%d\\n", unmapped, low, no_as, other > counts }' |
            samtools fastq -@ "${task.cpus}" -N - |
            tee fastq.fifo |
            bgzip -@ "${task.cpus}" -c > "\$READS"; then
            echo "ERROR: Could not decode ${alignment}. For CRAM, provide the matching --cram_reference if required." >&2
            exit 1
        fi
        wait \$P_PRIMARY
        wait \$P_FASTQ
        read -r UNMAPPED_RECORDS LOW_AS_READS NO_AS_RECORDS OTHER_RECORDS < low_as.classes
        if [[ "\$OTHER_RECORDS" != 0 || \$((UNMAPPED_RECORDS + LOW_AS_READS)) != "\$(< fastq.count)" ]]; then
            echo "ERROR: ${sampleID}: extraction classes do not add up: FLAG 4 \$UNMAPPED_RECORDS + low-AS \$LOW_AS_READS != \$(< fastq.count) FASTQ reads (unclassified: \$OTHER_RECORDS)" >&2
            exit 1
        fi
        touch low_as.ids
        gzip -c low_as.ids > "${sampleID}.low_as_read_ids.txt.gz"
        echo "extractReads: ${sampleID}: \$UNMAPPED_RECORDS FLAG 4 + \$LOW_AS_READS low-AS mapped reads (AS < 30); \$NO_AS_RECORDS mapped primary records had no AS tag and were not added" >&2
        rm -f fastq.fifo fastq.count low_as.classes low_as.ids
    else
        mkfifo unmapped.fifo
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
        UNMAPPED_RECORDS="\$(< unmapped.count)"
        rm -f unmapped.fifo unmapped.count
    fi
    # QC floor (pksProfiler COHORT_CHECKS A3; review #34). An input whose unmapped reads were removed
    # before it reached us -- the @PG chain shows samtools view -F 4 / -f 2 / -q / a region, or a
    # sambamba filter (input_unmapped_filter, scripts/header_provenance.py) -- or that has no FLAG 4
    # record at all cannot be profiled: its empty or rescue-only profile would read as "no microbes".
    # A real WGS/WES/RNA alignment always leaves some reads unmapped. Recorded here; main.nf acts on it.
    UPSTREAM_FILTER="\$(awk -F '\\t' 'NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "input_unmapped_filter") c = i } NR == 2 && c { print \$c }' "${sampleID}.input_provenance.tsv")"
    INPUT_QC_FLOOR=pass
    if [[ -n "\$UPSTREAM_FILTER" && "\$UPSTREAM_FILTER" != "none" ]]; then
        INPUT_QC_FLOOR=unmapped_filtered_upstream
    elif [[ "\$UNMAPPED_RECORDS" == 0 ]]; then
        INPUT_QC_FLOOR=no_unmapped_records
    fi
    if [[ "\$INPUT_QC_FLOOR" != pass ]]; then
        echo "extractReads: ${sampleID}: NOT ASSESSABLE (\$INPUT_QC_FLOOR; upstream filter: \${UPSTREAM_FILTER:-none}, FLAG 4 records: \$UNMAPPED_RECORDS)" >&2
    fi
    printf 'sample\tlibrary_primary_records\textracted_unmapped_records\tlow_as_reads_extracted\tlow_as_no_as_records\tinput_qc_floor\n%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${sampleID}" "\$(< primary.count)" "\$UNMAPPED_RECORDS" "\$LOW_AS_READS" "\$NO_AS_RECORDS" "\$INPUT_QC_FLOOR" > "${sampleID}.extract_counts.tsv"
    rm -f primary.fifo primary.count
    gzip -t "\$READS"
    """
}
