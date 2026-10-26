# What the host-depletion indexes contain

`--hg38_db` and `--t2t_phix_db` are the two minimap2 indexes every run depletes against before
anything is classified, so every CMPipeline count is what these two did not remove. On TSCC both are
the shared files under `/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/`, the same two
files pksProfiler uses. Both are rows in [`conf/databases.tsv`](../conf/databases.tsv), and every run
records whether they are still the frozen files (*DATABASE VERSIONS* in `RUN_REPORT.txt`).

The contents below were read out of the index files themselves (their header and sequence dictionary)
by pksProfiler on 2026-10-22 (`scripts/describe_minimap2_index.py`, documented in pksProfiler's
`docs/host_depletion_references.md`). They were not re-read for this page; `conf/databases.tsv`
carries pksProfiler's dictionary sha256 for each index as its key, so a run against a different index
reports `CHANGED: key`.

## `human-GRC-db.mmi` — `--hg38_db`

GRCh38.p14, the RefSeq full assembly (GCF_000001405.40): 705 sequences, 3,298,430,636 bp, built
`-x sr` (k=21, w=11). Primary chromosomes, mitochondrion, and the unlocalized, unplaced, ALT and patch
scaffolds. **Not** in it: hs38d1 decoys, chrEBV, HLA allele contigs, phiX.

## `human-GCA-phix-db.mmi` — `--t2t_phix_db`: there is no phiX in it

T2T-CHM13v2.0 (GCA_009914755.4): 25 sequences, 3,117,292,070 bp — the published CHM13v2.0 total to
the base, chr1–22, X, Y (HG002) and the mitochondrion. **phiX174 (NC_001422.1, 5,386 bp) is not one of
the 25 sequences**, despite the file name, the parameter name `--t2t_phix_db` and the QC columns
`reads_after_t2t_phix` / `pct_removed_t2t_phix`.

The code shows how that happened. `scripts/download_references.sh` (copied from
cguccione/human_host_filtration, Feb 2025) builds the index input with

```
cat ref/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna ref/NC_001422.fna > ref/human-GCA-phix.fna
```

but nothing in the script downloads `ref/NC_001422.fna`. `cat` writes the T2T FASTA, reports the
missing phiX file, and the script — which has no `set -e` — carries on; `create_minimap2_indexes.sh`
then indexes a T2T-only FASTA. That is consistent with the measured contents (T2T alone, to the base).
Whether the TSCC index was built by exactly this script is not recorded (no build log survives), so
this is the likely mechanism, not a proven history.

What it means:

- **What the pass does:** a second human-depletion pass against T2T-CHM13v2.0, removing reads from
  sequence GRCh38 does not represent. That is its main purpose and it is unaffected.
- **What it does not do:** remove phiX spike-in. phiX reads reach KrakenUniq and Bracken and **count in
  their totals and library sizes** (Bracken tables, `bracken_library_sizes.tsv`, decontamination library
  sizes). Viruses are dropped at the consensus step (`compute_consensus_taxa.py`), so consensus,
  decontaminated and batch-corrected tables contain no phiX taxon — unless `--skip_consensus`
  sends the Bracken tables straight on, in which case they do.
- **How to describe it:** in a methods section, "depleted against GRCh38.p14 and T2T-CHM13v2.0", not
  phiX depletion. Read `reads_after_t2t_phix` as "reads surviving the T2T pass".

**Resolution for now: documented, not changed.** The index, the parameter `--t2t_phix_db` and the
column names stay, because rebuilding the index mid-cohort would change host depletion for every
sample already run (port review #19; pksProfiler made the same call). The cohort report and docs now
label the pass "T2T-CHM13". Adding phiX (to the index, or as a tiny third pass) is a behaviour change
for a planned re-run boundary, with an A/B.

## Re-checking

Identity, cheap (a few kB of each index; what `conf/databases.tsv` freezes as `mmi:dictionary`): every
run does it. Integrity, expensive (reads the whole 9 GB file; run deliberately, not per run):

```
sha256sum /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi
sha256sum /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi
```

---
Back to the [README](../README.md).
