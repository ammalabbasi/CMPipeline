# batch_scripts

## Safe to submit

- `main.sbatch` — launches the pipeline with `-profile tscc`
- `resume_main.sbatch` — the same, kept under the name older notes use

        cd <your CMPipeline checkout>
        mkdir -p logs
        sbatch scripts/batch_scripts/main.sbatch samples.csv

  Both take the sample sheet as `$1` and pass any further arguments to Nextflow. They
  run in `$SLURM_SUBMIT_DIR`, so they operate on the checkout you submit from.

## Not safe to submit as they are

The eleven other `*.sbatch` files here **have no `#SBATCH` directives at all** — they
were stripped in commit `ca7b846` — and each one `cd`s into a hard-coded path under
another user's directory.

Two consequences, both silent:

1. No directives means Slurm applies the **1 hour default walltime**, so a database
   download or index build is killed partway with `TIMEOUT` and a partial file on disk.
2. The hard-coded `cd` means the job does its work in someone else's tree, or fails
   at the `cd` if you cannot read it.

| script | `#SBATCH` directives | hard-coded foreign path |
|---|---|---|
| build_metaphlan_index.sbatch | 0 | yes |
| create_minimap2_indexes.sbatch | 0 | yes |
| download_gtdb_checkm2.sbatch | 0 | yes |
| download_microbial_db.sbatch | 0 | yes |
| download_references.sbatch | 0 | yes |
| repair_pairs.sbatch | 0 | yes |
| repair_pairs_large.sbatch | 0 | yes |
| run_antismash.sbatch | 0 | yes |
| run_consensus_taxa.sbatch | 0 | yes |
| run_mag.sbatch | 0 | yes |
| unpack_microbial_dbs.sbatch | 0 | yes |

These are one-time setup jobs, so they have not been rewritten: the right walltime and
memory for each needs measuring, not guessing. Until then, run them by hand with explicit
`sbatch` arguments, e.g.

    sbatch -A ddp302 -p platinum -q hcp-ddp302 -c 4 --mem 32G -t 24:00:00 \
        --chdir "$PWD" scripts/batch_scripts/download_references.sbatch

and check where the script `cd`s before you do. Directives added to these files must sit
directly under the shebang — Slurm stops reading them at the first executable line.
Validate with `sbatch --test-only`.
