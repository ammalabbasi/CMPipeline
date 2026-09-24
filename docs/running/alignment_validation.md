# Alignment validation (opt-in)

`--run_alignment_validation true` adds a step between consensus and decontamination that asks
whether each taxon's reads actually come from its genome. It catches misassignment that two
agreeing classifiers can share: reads piled on rRNA operons or a plasmid, reads of an absent
relative, and "shadow" species that only share reads with a real one.

**When to turn it on.** For a cohort whose species-level calls you will report or build on,
once the cohort has run through consensus at least once and you know which taxa matter. It is off
by default because it is the most expensive step after host depletion: it downloads and indexes
RefSeq genomes for every species in the table and re-aligns every sample's host-depleted reads.
It needs the full run (`--start_from beginning`, host depletion and classification on), the
Kraken `taxDB` (from `--kraken_db`, or `--validation_taxdb`), and internet on the compute node
unless `--validation_ref` supplies a pre-built reference.

**What changes downstream.** Bracken counts remain the quantity. Validation only sets a count to
0 where the taxon *failed* in that sample, and decontamination then runs twice, on the validated
species table and on its genus roll-up.

1. **Species list.** Species in the table going downstream (consensus, or Bracken when consensus
   passed through) with >= `--validation_min_support` reads in at least one sample, matched to NCBI
   taxids through the Kraken database's own `taxDB`.
2. **Reference**, built once per species list and cached under `--validation_cache_dir`:
   - up to `--validation_max_genomes` (50) RefSeq genomes per species (Complete Genome if any), de-duplicated at
     `--validation_dedup_ani` (99.5%) with skani
   - **human T2T in the same index**, so leftover human reads compete
   - one single-part minimap2 index (the run stops if it would split)
   - every accession and checksum in `validation_reference_manifest.tsv`
   - needs internet on the compute node
3. **Alignment.** Each sample's host-depleted reads, `--validation_batch_size` samples per task in one
   minimap2 stream, so the index loads once per batch.
4. **Status** per sample x species, with every statistic in `06_VALIDATION/validation_stats.all.tsv`:
   - tested in this order: `unsupported` if Bracken gave it at least 25 reads but fewer than 10 align to
     it (`--validation_min_confirmed`; the cutoff PRISM uses; the reads are not this species');
     `insufficient` if fewer than 25 supporting reads over all its genomes;
     `failed` if the unique fraction is `< 0.02` (its reads all belong to another species too);
     `insufficient` if fewer than 25 reads land on its best genome (spread thin over many strains);
     `validated` if the 10 kb spread `bins_ratio >= 0.5` and median identity `>= 0.98`; otherwise `failed`,
     with `fail_reason` naming each test
   - `not_in_reference`: no genome in the reference, or Homo_sapiens (never scored)
   - pooled over all of a species' genomes: supporting/unique/multi reads, unique fraction, identity;
     best genome only: breadth, bins, `reads_on_best`. `minimap2_count` equals `unique_reads`.
   - E. coli + Shigella are scored as one taxon (`ref/validation_groups.tsv`).
5. **Mask.** Bracken counts stay the quantity. A count is set to 0 **only where the taxon failed or was
   unsupported** in that sample. The masked species table and its genus roll-up both go to decontamination, and
   `decontam_species_vs_genus.tsv` compares the two sets of calls. Batch correction continues from
   `--taxonomic_rank`.

The thresholds come from a simulation benchmark on public genomes (21 CRC-relevant species; reliable from
about 100 read pairs per taxon, zero false calls). **Identity 0.98 is the least certain for real reads**
(about 1% error); check `validation_stats.all.tsv` after the first real run. `--validation_export_reads true`
also writes the reads supporting each validated taxon. Those are sample-derived reads: protected data.

| Parameter | Default | Meaning |
|---|---|---|
| `--run_alignment_validation` | `false` | turn the lane on |
| `--validation_min_support` | `25` | reads on the best genome for a species to be listed, and for a sample x species to be scored at all |
| `--validation_min_bins_ratio` | `0.5` | 10 kb windows hit / windows expected for that many read pairs |
| `--validation_min_identity` | `0.98` | median alignment identity |
| `--validation_min_unique_frac` | `0.02` | reads unique to the species / all its reads |
| `--validation_max_genomes` | `50` | RefSeq genomes per species, before de-duplication |
| `--validation_dedup_ani` | `99.5` | skani ANI at which two genomes of a species count as one |
| `--validation_max_missing_frac` | `0.10` | the build fails if more than this fraction of species get no genome |
| `--validation_batch_size` | `20` | samples per alignment task |
| `--validation_groups` | `ref/validation_groups.tsv` | species scored together as one taxon (E. coli + Shigella) |
| `--validation_taxdb` | `<kraken_db>/taxDB` | the NCBI taxonomy used to match names to taxids |
| `--validation_cache_dir` | `.validation_cache/` in the repo, or `$CMP_VALIDATION_CACHE` | where built references are kept and reused |
| `--validation_ref` | none | a pre-built `ref-<key>/` folder from the cache: skips the species list and the build (no internet needed; pins one reference across cohorts) |
| `--validation_min_confirmed` | 10 | a species with at least `--validation_min_support` Bracken reads but fewer supporting alignments than this is `unsupported`, and its count is zeroed (decision 29; PRISM's cutoff) |
| `--validation_max_ref_gb` | 40 | refuse a planned reference larger than this many Gb of sequence, checked from RefSeq genome sizes **before** anything is downloaded or indexed |
| `--validation_refseq_snapshot` | the month the run starts | RefSeq snapshot month (`YYYY-MM`), part of the reference's cache key. Within a month the cached reference is reused; a new month builds a fresh, correctly labelled one. A past month is reused only if still cached; otherwise the build stops (pin it with `--validation_ref`) |
| `--validation_export_reads` | `false` | also write each validated taxon's supporting reads (protected data) |

Fractions must be in [0, 1] and counts positive integers; both are checked at launch.

**Outputs**, under `06_VALIDATION/` (`--validation_dir`): `validation_species.tsv` (the species
list), `validation_stats.all.tsv` (every statistic per sample x species), `validation_matrix.tsv`
(the status grid), `validated.species.tsv` and `validated.genus.tsv` (the masked tables that go to
decontamination), `validation_mask_report.tsv` (per species: samples validated, failed and insufficient,
and reads removed; also shown in the cohort report), `decontam_species_vs_genus.tsv`, and `per_sample/` with each batch's statistics. The
reference itself, with `validation_reference_manifest.tsv` and `missing_species.tsv`, is in the
cache folder.

---
[Back to the README](../../README.md)
