nextflow.enable.dsl=2

process filterReads {
    label 'per_sample'   // --sample_failure_strategy (audit C10)
    tag "$sampleID"
    scratch true
    label 'filter_reads'
    // Intermediate read/alignment files are published only with --save_intermediates true (audit C11);
    publishDir("${params.unmapped_bam_dir}", mode: 'copy', saveAs: { fn -> (params.save_intermediates.toString() == 'true' || !(fn =~ /\.(fastq|fq|fasta|fa|sam|bam|bowtie2\.bz2)(\.gz|\.bz2)?$/)) ? fn : null })
    conda "${params.fastp_env}"
    maxRetries 2

    input:
    tuple val(sampleID), path(fastq_files), val(trim_adapters)

    // Named emits: a second output turns the process handle into a multi-channel object,
    // and `filterReads(...).set{}` at the call site would fail with "Multi-channel output
    // cannot be applied to operator set".
    // Digests of the helper scripts / references this task uses (Provenance.code/.data). A val
    // input, so editing a script or swapping a reference invalidates -resume for exactly the
    // dependent tasks; embedding it in the script text does NOT (verified 2026-09-23, audit C03).
    val deps
    output:
    tuple val(sampleID),
        path("${sampleID}.UNMAPPED.FASTP.FILTERED.fastq.gz"),
        path("${sampleID}.reads_after_filtering.txt"), emit: reads
    // fastp writes these either way; before this they were produced and then discarded,
    // throwing away the richest QC signal in the pipeline. MultiQC reads the JSON.
    path("${sampleID}.fastp.json"), emit: fastp_json
    path("${sampleID}.fastp.html"), emit: fastp_html

    script:
    def input_list = fastq_files instanceof Collection ? fastq_files : [fastq_files]
    def inputs = input_list.collect { "\"${it}\"" }.join(' ')
    def tagged_commands = input_list.withIndex().collect { input_fastq, index ->
        def mate = index + 1
        """gzip -dc "${input_fastq}" | awk -v mate="${mate}" '
        NR % 4 == 1 {
            split(\$0, parts, " "); id=parts[1]
            if (id !~ /\\/[12]\$/) id=id "/" mate
            \$0=id substr(\$0, length(parts[1]) + 1)
        }
        { print }'"""
    }.join('\n')
    def fastp_input = input_list.size() == 1 ? "-i ${inputs}" : "--stdin"
    def fastp_prefix = input_list.size() == 1 ? "" : "{\n${tagged_commands}\n} |"
    // Adapter scanning is decided per sample in main.nf (see --adapter_trim).
    def adapter_args = trim_adapters ? "--adapter_fasta \"${params.adapters}\"" : ""

    """
    set -euo pipefail
    FILTERED="${sampleID}.UNMAPPED.FASTP.FILTERED.fastq.gz"
    COUNT="${sampleID}.reads_after_filtering.txt"
    for input_fastq in ${inputs}; do
        gzip -t "\$input_fastq"
    done
    ${fastp_prefix} fastp -l 45 ${adapter_args} --cut_tail \
        ${fastp_input} -w "${task.cpus}" --compression 4 \
        --json "${sampleID}.fastp.json" --html "${sampleID}.fastp.html" \
        -o "\$FILTERED"

    # fastp exits 0 on an input with no reads but leaves a 0-byte output, which is not
    # valid gzip. Write a real empty member so the file downstream is always readable;
    # the read count below is what decides whether the sample continues.
    if [[ ! -s "\$FILTERED" ]]; then
        : | gzip -c > "\$FILTERED"
    fi
    gzip -t "\$FILTERED"

    awk '
        /"after_filtering"[[:space:]]*:/ { section=1; next }
        section == 1 && /"total_reads"[[:space:]]*:/ {
            value=\$0
            gsub(/[^0-9]/, "", value)
            print value
            exit
        }
    ' "${sampleID}.fastp.json" > "\$COUNT"
    [[ -s "\$COUNT" ]] || echo 0 > "\$COUNT"
    """
}
