# June-Young integration v2 (checkpoint 3, not a packaged release)

Baseline: c39a8f43ed6661679e4254abdd712a6ca832600b. The v1 directory is preserved.

- Bracken preprocessing emits canonical rank-specific tables and a separate validated library-size table. Raw hierarchical MPA is not an accepted library-size substitute.
- Classification can start from the pipeline's host-depleted FASTQ with --skip_host_depletion true and a patient,fastq_1 CSV. One normalized compressed FASTQ per sample is required. Host mapping/filtering are not repeated.
- In full or host-depleted-read runs, decontamination acts on species. Genus counts are then summed from surviving species, without a second genus test; a species with no genus parent stays in the species table and is listed in <prefix>.genus_rollup_excluded.tsv instead of stopping the run. Batch correction runs on both the species table and the derived genus table (05_BATCH_CORRECTION/<run>_species and <run>_genus). Original species and derived genus decontamination outputs are both published under species_run.
- Standalone --start_from decontam treats --taxonomic_rank as the input rank. A species input also emits derived genus output; an explicitly supplied genus input remains a separate standalone genus analysis, not an equivalent species-first run.
- Library-size values and alignment statuses are strict inputs. Unknown/malformed results are errors; legitimate insufficient/reference-missing alignment states keep the existing policy.
- Retained uncallable taxa are marked retained_not_callable. An all-uncallable analysis reports insufficient_evidence, not a scientific clearance. The existing skip/error policy is unchanged.
- All-zero policy, alignment thresholds and existing consensus pass-through policy are unchanged.
- Run tests/integration/run_checkpoint3.sbatch on a compute node. It checks bounded Python/R behaviors, executes the official species decontamination/aggregation path, and previews the host-depleted entry. It does NOT execute real upstream classification, build fresh environments, perform CRC12 validation, or validate batch correction.
- The general README will be reconciled at the packaging checkpoint. This working copy is not yet ready for release.
