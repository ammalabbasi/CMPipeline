# Decontamination

Batch-aware, prevalence-based contamination detection with the `decontam` R package
(`scripts/decontamination.R`). It is reached three ways: as part of a full run
(`--run_decontam true`), from an existing consensus table
(`--start_from decontam --consensus_otu_table <path>`), or not at all (`--run_decontam false`).
**`--run_decontam` defaults to `true`**, so a full run must either supply the parameters below
and `--metadata_file`, or pass `--run_decontam false`; otherwise it stops at launch.

**Four parameters are required whenever decontamination runs, and have no default.**
Naming which samples are the real material and which stand in for blanks is a scientific
choice, not something the pipeline should guess:

| Parameter | Meaning |
|---|---|
| `--taxonomic_rank` | `genus` or `species` — which consensus/Bracken table to decontaminate |
| `--control_role` | `experimental_blank` (true blanks) or `surrogate_biological` (a sample class standing in for blanks) |
| `--positive_values` | the sample types that are real material, i.e. **not** controls |
| `--control_values` | the sample types used as negative controls |

A run that omits any of them stops in about a second with a named error rather than
guessing. `--type_column` says which metadata column holds the role and defaults to
`sample_type`.

For TCGA, one flag supplies the whole contract — the manuscript's model, with blood-derived
normals as surrogate negative controls:

```bash
nextflow run main.nf -profile tscc --sample samples.csv \
  --metadata_file metadata.tsv --run_decontam true \
  --cohort_preset tcga --taxonomic_rank genus
```

`--cohort_preset tcga` sets `--type_column "Sample Type_x"`,
`--positive_values "Primary Tumor,Solid Tissue Normal"` and
`--control_values "Blood Derived Normal"`. Passing any of them a conflicting value is an
error rather than a silent override. Note that the positives include Solid Tissue Normal:
it is real material, not a blank.

For a cohort with Tumor/Normal in a `Type` column, state them:

```bash
nextflow run main.nf -profile tscc --sample samples.csv \
  --metadata_file metadata.tsv --run_decontam true --taxonomic_rank genus \
  --type_column Type --control_role experimental_blank \
  --positive_values Tumor --control_values Normal
```

A sample whose type matches neither list stops the run, because
`--unselected_type_policy` defaults to `error`; set it to `exclude` to drop those samples
instead.

The other columns and table layout:

| Parameter | Default | Meaning |
|---|---|---|
| `--metadata_file` | none (required) | per-sample metadata, tab- or comma-separated |
| `--sample_id_column` | `sampleid` | metadata column holding the sample IDs |
| `--batch_column` | `shipment_batch` | metadata column holding the batch; also used by batch correction |
| `--taxon_column` | `clade_name` | taxonomy column in the OTU table; every other column is a sample |

The ID, type and batch columns are checked against the metadata header at launch.

**How a contaminant is called.** Taxa are pre-filtered, then tested batch by batch with
decontam's prevalence method, then called across batches:

| Parameter | Default | Meaning |
|---|---|---|
| `--decontam_min_abundance` | `5` | a taxon's reads in a sample below this count as absent |
| `--decontam_min_prevalence` | `0.05` | after that filter, a taxon must be present in at least 5% of samples to be tested |
| `--decontam_threshold` | `0.1` | decontam's prevalence-test threshold |
| `--decontam_min_batches` | `2` | a taxon flagged in at least this many batches is a contaminant |
| `--min_total_per_batch`, `--min_positive_per_batch`, `--min_control_per_batch` | `5`, `1`, `1` | a batch with fewer samples, positives or controls than this is not tested |

`--decontam_min_prevalence` 0.05 matches the published method ("≥5 reads, detected in ≥5% of
samples" for CRC). It was 0.02 until 2026-09-22; results produced before then are not
comparable.

Contaminants are removed **globally**, not zeroed per batch.

**All-zero samples.** A sample whose counts are all zero in the table being decontaminated
(for example, a normal that carries none of the consensus genera) is evidence that those taxa
are absent, not a broken sample:

| Parameter | Default | Meaning |
|---|---|---|
| `--decontam_zero_total_policy` | `keep` | `keep`: an all-zero sample stays in the prevalence test, as long as its library was not empty. `drop`: every all-zero sample is excluded (the older behaviour, which could remove a batch's controls) |
| `--decontam_min_library` | `1000` | under `keep`, an all-zero sample whose total Bracken library is below this is excluded as `insufficient_library` |
| `--decontam_library_sizes` | none | the library-size table for `--start_from decontam`: pass the `bracken_library_sizes.tsv` that consensus wrote. In a full run it is supplied automatically |
| `--decontam_allow_unknown_library` | `false` | an all-zero sample whose library size is unknown stops the run unless this is `true`, in which case it is kept |

**Taxa that cannot be tested.** A taxon present in a batch but in fewer than 2 of its samples
cannot be tested there; decontam returns no answer.

| Parameter | Default | Meaning |
|---|---|---|
| `--decontam_unevaluable_policy` | **`skip`** | `skip`: record the taxon as `unevaluable` in that batch and not flagged there; the `--decontam_min_batches` rule then decides from the batches that could test it. `error`: stop the run |

Under `skip`, a taxon testable in fewer than `--decontam_min_batches` batches can never be called
a contaminant, however contaminant-like it looks. `contaminants_final.tsv` shows this per taxon
in `n_evaluable_batches` and `callable`; read "not a contaminant" as "could not be called" where
`callable` is `FALSE`.

**Outputs**, under `04_DECONTAMINATION/<run>/`, all prefixed with the run name:
`decontam_pkg_decontaminated.csv` (the cleaned table), `contaminants_by_batch.tsv`,
`contaminants_final.tsv`, a per-taxon `filter_report.tsv` carrying a `terminal_reason`,
`sample_validation.tsv` (every sample's role and whether it was used or excluded, and why),
`validated_metadata.tsv`, `batch_eligibility.tsv`, `run_manifest.tsv` (every setting and input
checksum), `decontamination_summary.txt` and `decontamination_plots/`.

To decontaminate an existing consensus table without rerunning anything upstream:

```bash
nextflow run main.nf -profile tscc --start_from decontam --run_decontam true \
  --consensus_otu_table CONSENSUS_TAXA/bracken.metaphlan.common.genus.mpa.report.txt \
  --decontam_library_sizes CONSENSUS_TAXA/bracken_library_sizes.tsv \
  --metadata_file metadata.tsv --taxonomic_rank genus --cohort_preset tcga
```

---
[Back to the README](../../README.md)
