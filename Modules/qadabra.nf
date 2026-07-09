nextflow.enable.dsl=2

process Qadabra {

    tag "${prefix}_${level}_${model}"
    label 'process_high'

    conda "${params.qadabra_env}"
    publishDir "${params.qadabra_dir}/${prefix}_${level}_${model}", mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata), val(level), val(model), val(formula)

    output:
    tuple val("${level}_${model}"), path("results/*/concatenated_differentials.tsv"), optional: true, emit: results
    tuple val("${level}_${model}"), path("results/*/tools/*/differentials.tsv"),       optional: true, emit: per_tool

    script:
    def dataset_name = "${prefix}_${level}"
    def factor_col   = params.type_column ?: "Type"
    """
    # Convert TSV count table to BIOM (required by Qadabra)
    python - <<'PYEOF'
import pandas as pd
from biom.table import Table as BiomTable
from biom.util import biom_open

df = pd.read_csv("${otu_table}", sep="\\t", index_col=0)
df.index = df.index.astype(str)

# Filter to level
if "${level}" == "genus":
    cols = [c for c in df.columns if c.split("|")[-1].startswith("g__")]
else:
    cols = [c for c in df.columns if "|s__" in c]
df = df[cols]

# Remove zero taxa and round to int
df = df.loc[:, df.sum() > 0]
mat = df.T.values.astype(int)
table = BiomTable(mat, observation_ids=df.columns.tolist(), sample_ids=df.index.tolist())
with biom_open("table.biom", "w") as f:
    table.to_hdf5(f, "qadabra")
print(f"Wrote BIOM: {mat.shape[0]} taxa x {mat.shape[1]} samples")
PYEOF

    # Create Qadabra workflow directory
    qadabra workflow create \\
        --table table.biom \\
        --metadata ${metadata} \\
        --factor-name ${factor_col} \\
        --reference-factor-level ${params.da_reference} \\
        --dataset-name ${dataset_name} \\
        --verbose

    # Run Snakemake (conda envs pre-built in the Qadabra conda prefix)
    snakemake \\
        --use-conda \\
        --cores ${task.cpus} \\
        --rerun-incomplete \\
        --keep-going \\
        --printshellcmds \\
        2>&1
    """
}
