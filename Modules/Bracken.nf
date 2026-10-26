nextflow.enable.dsl=2

process Bracken {
    label 'per_sample'   // --sample_failure_strategy (audit C10)

    scratch true
    label 'process_high_disk'
    // Intermediate read/alignment files are published only with --save_intermediates true (audit C11);
    // metaphlan.fastq.gz (a copy of the reads) and the classified/unclassified FASTA included.
    publishDir "${params.krakenuniq_bracken_dir}", mode: 'copy', saveAs: { fn -> (params.save_intermediates.toString() == 'true' || !(fn =~ /\.(fastq|fq|fasta|fa|sam|bam|bowtie2\.bz2)(\.gz|\.bz2)?$/)) ? fn : null }
    conda "${params.krakenuniq_bracken_env}"

    input:
    tuple val(sampleID), path(reads)

    // `optional: true` does NOT work inside a tuple: Nextflow still tries to resolve
    // every component and fails with "ls: cannot access '*.metaphlan.fastq.gz'" when the
    // glob matches nothing. That is precisely the designed-for case -- a sample below
    // --consensus_min_reads has no MetaPhlAn input -- so the gate killed the run whenever
    // it actually fired. The optional file is therefore its own named output.
    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple val(sampleID), path("*.krakenuniq.report.txt"),
      path("*.classified.fasta"),
      path("*.unclassified.fasta"),
      path("*.bracken.*.report.txt"),
      path("*.bracken.*.krakenreport.txt"),
      path("*.bracken.*.mpa.krakenreport.txt"), emit: reports
    tuple val(sampleID), path("*.bracken_read_length.tsv"), emit: read_length
    tuple val(sampleID), path("*.metaphlan.fastq.gz"), emit: metaphlan_reads, optional: true

    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = sampleID
    def read_files = reads instanceof Collection ? reads : [reads]
    def inputs = read_files.collect { "\"${it}\"" }.join(' ')

    """

    KRAKEN_REPORT="${SAMPLE_NAME}.krakenuniq.report.txt"
    CLASSIFIED="${SAMPLE_NAME}.classified.fasta"
    UNCLASSIFIED="${SAMPLE_NAME}.unclassified.fasta"
    OUTPUT="${SAMPLE_NAME}.krakenuniq.output.txt"
    BOUT_G="${SAMPLE_NAME}.bracken.G.report.txt"
    BOUT_S="${SAMPLE_NAME}.bracken.S.report.txt"

    # Actual execution

    echo "SAMPLE_NAME is: $SAMPLE_NAME"
    gzip -dc ${inputs} > "${SAMPLE_NAME}.trimmed.fastq"

    CLASSIFIED_FASTA="${SAMPLE_NAME}.classified.fasta"
    UNCLASSIFIED_FASTA="${SAMPLE_NAME}.unclassified.fasta"
    REPORT="${SAMPLE_NAME}.krakenuniq.report.txt"
    OUTPUT="${SAMPLE_NAME}.krakenuniq.output.txt"

    krakenuniq --db "${params.kraken_db}" --threads ${task.cpus} --report-file \${REPORT} --output \${OUTPUT} \
        --classified-out \${CLASSIFIED_FASTA} --unclassified-out \${UNCLASSIFIED_FASTA} "${SAMPLE_NAME}.trimmed.fastq"

    # Bracken read length (2026-10-26; pksProfiler c6599c8 measures it per sample too). Bracken's
    # re-estimation uses a kmer distribution built for one read length; a fixed 50 was applied to
    # ~150 bp TCGA WGS reads. auto: median length of the first 1,000,000 reads, then the database
    # length closest to it (ties -> shorter). The choice is written out for QC and provenance.
    MEDIAN_LEN=\$(awk 'NR%4==2{print length(\$0)} NR>=4000000{exit}' "${SAMPLE_NAME}.trimmed.fastq" | sort -n | awk '{a[NR]=\$1} END{print (NR ? a[int((NR+1)/2)] : 0)}')
    AVAILABLE=\$(ls "${params.kraken_db}" | sed -n 's/^database\([0-9][0-9]*\)mers\.kmer_distrib\$/\1/p' | sort -n | tr '\n' ' ')
    if [[ "${params.bracken_read_length}" == "auto" ]]; then
        BRACKEN_LEN=\$(for L in \$AVAILABLE; do d=\$(( L > MEDIAN_LEN ? L - MEDIAN_LEN : MEDIAN_LEN - L )); echo "\$d \$L"; done | sort -n -k1,1 -k2,2 | head -1 | cut -d' ' -f2)
        LEN_MODE=auto
    else
        BRACKEN_LEN=${params.bracken_read_length}
        LEN_MODE=fixed
    fi
    [[ -n "\$BRACKEN_LEN" ]] || { echo "ERROR: no database<N>mers.kmer_distrib in ${params.kraken_db}" >&2; exit 1; }
    printf 'sample\tmedian_read_length\tbracken_read_length\tmode\tavailable_lengths\n%s\t%s\t%s\t%s\t%s\n' \
        "${SAMPLE_NAME}" "\$MEDIAN_LEN" "\$BRACKEN_LEN" "\$LEN_MODE" "\$(echo \$AVAILABLE | tr ' ' ',')" > "${SAMPLE_NAME}.bracken_read_length.tsv"

    for lvl in G S; do
        bracken_output="${SAMPLE_NAME}.bracken.\${lvl}.report.txt"
        bracken_kraken_report="${SAMPLE_NAME}.bracken.\${lvl}.krakenreport.txt"
        bracken_kraken_mpa_report="${SAMPLE_NAME}.bracken.\${lvl}.mpa.krakenreport.txt"
        # A thin sample where no taxon at this level reaches -t makes Bracken print "Error: no reads
        # found" and exit 1, which under errorStrategy 'finish' ended the whole cohort (audit C02).
        # That one case becomes an explicit below-threshold no-call with the normal file layout
        # (header-only report and mpa); any other Bracken error still fails the task.
        set +e
        bracken -d "${params.kraken_db}" -i \${REPORT} -o \${bracken_output} -w \${bracken_kraken_report} -r \${BRACKEN_LEN} -l \${lvl} -t ${params.bracken_threshold} 2> bracken.\${lvl}.err
        bracken_rc=\$?
        set -e
        cat bracken.\${lvl}.err >&2
        if [[ \$bracken_rc -ne 0 ]]; then
            if grep -q "no reads found" bracken.\${lvl}.err; then
                echo "Bracken \${lvl}: below_threshold for ${SAMPLE_NAME} (no taxon with >= ${params.bracken_threshold} reads); writing an empty no-call"
                printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n' > \${bracken_output}
                : > \${bracken_kraken_report}
                printf '#Classification\t%s\n' "\$(basename \${bracken_kraken_report})" > \${bracken_kraken_mpa_report}
                continue
            fi
            exit \$bracken_rc
        fi
        python "${params.scripts}"/kreport2mpa.py -r \${bracken_kraken_report} -o \${bracken_kraken_mpa_report} --display-header
    done

    # Gate MetaPhlAn per sample using the classified microbial read count.
    #
    # This read three things wrong and silently returned 0 for every sample from
    # 2026-09-04 until 2026-09-23, so MetaPhlAn never ran and consensus never had
    # MetaPhlAn data:
    #   1. it took the rank from \$4, but a KrakenUniq report has 9 columns
    #      (% reads taxReads kmers dup cov taxID rank taxName) -- \$4 is kmers;
    #   2. the rank of both rows is the literal string "no rank", never "R" or "U",
    #      so no rank comparison can ever match. They are identified by taxName;
    #   3. root and unclassified are SIBLINGS, so root-minus-unclassified is the wrong
    #      arithmetic. Root's own clade count is the classified total.
    # taxName is trimmed before comparing, exactly as build_qc_summary.py does (audit A16), so the
    # gate and the QC column cannot disagree on a padded name.
    MICROBIAL_READS=\$(awk -F '\t' '
        !/^#/ {
            n=\$9; gsub(/^[ \t\r]+|[ \t\r]+\$/, "", n)
            if(n=="root") {
                found++; v=\$2; gsub(/^[ \t\r]+|[ \t\r]+\$/, "", v)
                if(NF!=9 || v !~ /^[0-9]+\$/) bad=1
                total=v
            }
        }
        END {
            if(found!=1 || bad) { print "ERROR: expected exactly one valid KrakenUniq root count" > "/dev/stderr"; exit 1 }
            print total
        }' "\${REPORT}")
    echo "Bracken microbial reads for ${SAMPLE_NAME}: \${MICROBIAL_READS}"
    if [[ "\${MICROBIAL_READS}" -ge "${params.consensus_min_reads}" ]]; then
        gzip -c "${SAMPLE_NAME}.trimmed.fastq" > "${SAMPLE_NAME}.metaphlan.fastq.gz"
        echo "MetaPhlAn eligible: ${SAMPLE_NAME} (\${MICROBIAL_READS} >= ${params.consensus_min_reads})"
    else
        echo "MetaPhlAn skipped: ${SAMPLE_NAME} (\${MICROBIAL_READS} < ${params.consensus_min_reads})"
    fi

    echo "Done processing ${SAMPLE_NAME}"

    """
}
