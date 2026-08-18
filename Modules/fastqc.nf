nextflow.enable.dsl=2

process FASTQC {
    scratch true
    label 'fastqc'
    publishDir("${params.fastqc_dir}", mode: 'copy')
    conda "${params.fastqc_env}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(r1_fastq), path(r2_fastq)

    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:

    """
    # Defined report naming convention (FastQC: <basename>_fastqc.html/.zip)
    # Process basename using shell
    r1_base=\$(basename ${r1_fastq.toString()} .fastq.gz)
    r2_base=\$(basename ${r2_fastq.toString()} .fastq.gz)
    out1="\${r1_base}_fastqc.html"
    out2="\${r1_base}_fastqc.zip"
    out3="\${r2_base}_fastqc.html"
    out4="\${r2_base}_fastqc.zip"

    # Skip check
    if [[ -f "${params.fastqc_dir}/\$out1" && -f "${params.fastqc_dir}/\$out2" && -f "${params.fastqc_dir}/\$out3" && -f "${params.fastqc_dir}/\$out4" ]]; then
        echo "Skipping FASTQC: Found all existing fastqc outputs"

        for f in "\$out1" "\$out2" "\$out3" "\$out4"; do
            if [[ ! -f "\$f" ]]; then
                ln -s "${params.fastqc_dir}/\$f" . 2>/dev/null || cp "${params.fastqc_dir}/\$f" .
            fi
        done
        exit 0
    fi

    # Actual FastQC execution
    # Running FastQC on the extracted R1 and R2 fastq.gz files
    fastqc -t 8 -o ./  ${r1_fastq.toString()} ${r2_fastq.toString()}
    """
}