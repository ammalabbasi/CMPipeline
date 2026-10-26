# Classification and consensus

How KrakenUniq/Bracken, the MetaPhlAn gate and consensus behave, plus adapter trimming and empty samples.

## Adapter trimming

`--adapter_trim` chooses the adapter FASTA fastp scans reads against:

- `auto` (default): raw FASTQ against `--adapters` (`ref/known_adapters.fna`, 234
  sequences); reads extracted from a BAM or CRAM against `--adapters_core`
  (`ref/illumina_core_adapters.fna`: TruSeq read 1/2, TruSeq2 SE/PE, Nextera).
- `always`: the 234-sequence list for every sample. `never`: no FASTA.

Without a FASTA fastp still auto-detects one adapter per read, but no more.

**Why BAM/CRAM reads are scanned at all.** Until 2026-10-26 `auto` skipped the scan for
them, on the assumption that aligned reads were trimmed before alignment. GDC did not trim
TCGA before `bwa mem`, which soft-clips adapter read-through, so the extracted reads still
carry it. Untrimmed adapter is classified by KrakenUniq as *synthetic construct* (UniVec),
which dominated Bracken output in CRC50 for one sequencing centre.

**Why not the full list.** It is BBTools' generic `adapters.fa`: about 200 of its 234
entries are the same TruSeq, RNA and Nextera adapters repeated once per 6-8 bp index
barcode, plus retired kits (SOLiD, DpnII/NlaIII gene expression). Matching starts at the
shared prefix, so the barcode variants add almost nothing, and fastp's cost scales with the
list: on one aligned CRAM `filterReads` took 2 h 50 m with it versus 3 m 15 s without. On
300,000 synthetic 151 bp reads with 60,000 planted adapter tails (fastp 1.3.7): no FASTA 5 s
and 30,470 reads trimmed; core list 17 s, 61,426; full list 645 s, 66,112.

The QC summary's `adapter_list` column records which FASTA each sample used (read from fastp's own
command line; fastp 1.x writes an `adapter_cutting` block even for auto-detection alone).

## Low-score mapped reads (TCGA GDC BAMs)

For BAM/CRAM input the pipeline classifies the reads the upstream aligner did not place on the
human reference. Until now that meant FLAG 4 (unmapped) records only. `bwa mem` reports a read as
unmapped only when its best local alignment scores below `-T` (default 30). **GDC aligned TCGA with
`bwa mem -T 0`**, so a microbial read with a chance ~19 bp match to human comes out *mapped* (low
`AS`, mostly soft-clipped, MAPQ 0-9) and FLAG 4 extraction never sees it. pksProfiler measured the
effect for pks reads (synthetic): TCGA kept 57-66% of them, `-T 30` 98.5%.

`--extract_low_score_mapped` (default `auto`) also extracts **primary mapped records with
`AS:i` < 30**, the set a default `bwa mem` would have left unmapped:

- `auto` turns it on only when the input's `@PG` header runs `bwa mem` or `bwa-mem2 mem`
  (also inside a wrapper such as `sentieon bwa mem`) **with `-T` < 30**. STAR, minimap2, bowtie2,
  `bwa aln`/`sampe`, `bwa mem` at the default or a higher `-T`, an unreadable `-T`, no `@PG`, or only post-processing `@PG` (samtools, Picard)
  leave it off: other aligners' `AS` is not on bwa's scale, and with `-T` >= 30 there is nothing
  to rescue. The decision is made per sample from the header alone; the command line is read in
  memory, never written.
- `on` applies the rule to every BAM/CRAM; `off` restores FLAG 4 only. FASTQ input is unaffected.
- Secondary (0x100) and supplementary (0x800) records are never extracted; mapped records without
  an `AS` tag are never extracted (only counted); FLAG 4 records are extracted as before.
- Rescued reads go through the same `samtools fastq -N`, so names, `/1` `/2` suffixes and
  reverse-strand orientation are exactly as for unmapped reads. With the rule off the extraction
  pipeline is unchanged.

What each sample got is recorded in the QC summary (`upstream_aligner`, `upstream_aligner_version`,
`upstream_bwa_min_score`, `low_as_mode`, `low_as_mode_reason`) along with how many rescued reads
were extracted, survived host depletion and were classified (`low_as_reads_extracted`,
`low_as_reads_after_depletion`, `low_as_reads_classified`, `low_as_share_of_bracken`; see
[outputs](../outputs.md)). Survival is counted by read name: `extractReads` writes the rescued names
(an intermediate, published only with `--save_intermediates true`) and `lowAsSurvival` looks them up
in the host-depleted FASTQ and in KrakenUniq's classified FASTA.

**Known limits.** Bracken re-estimates abundances from the KrakenUniq report and does not attribute
individual reads, so the rescued share of each *taxon* is not measured; `low_as_share_of_bracken`
is the rescued share of Bracken's input. A rescued read that host depletion renames is not matched:
minimap2 (`-ax sr`) pairs adjacent records with the same base name and `samtools fastq` then adds
`/1` `/2`, which leaves `/1` `/2` names unchanged but would rename an unpaired, suffix-less read that
happens to sit next to a same-named one (an undercount, never an overcount of survivors).

**This changes results for TCGA GDC BAMs** (and any `bwa mem -T` < 30 input): more reads are
classified. Whether the recovery is the same for every taxon, and whether it adds as much to blood
normals as to tumours, is being tested before cohort numbers are compared (test plan,
2026-10-01). To reproduce an older run exactly, pass `--extract_low_score_mapped off`.

## Samples with no reads

A sample that has no reads left after filtering is logged and dropped; the rest
of the run continues. This covers aligned input that carries no unmapped records,
which is common for RNA-seq BAMs written without STAR's `--outSAMunmapped Within`.

## Classification and consensus

KrakenUniq and Bracken run first for every sample. The sample's classified microbial read
count, the clade count of the `root` row of its KrakenUniq report, is a per-sample gate:
**MetaPhlAn 4 runs only when that count is at least `--consensus_min_reads` (100000).**
Samples below the gate are logged and skip MetaPhlAn, and therefore HUMAnN3 too (see
*Functional profiling*). Their Bracken counts still go forward. The QC summary's
`kraken_classified_reads` and `cleared_consensus_gate` columns show, per sample, the count and
whether it cleared. That count includes reads classified as human (`kraken_human_reads`, taxID
9606) and as synthetic construct (`kraken_synthetic_reads`, 32630: vectors, adapters, primers);
`bracken_microbial_reads` is the classified count minus both. `bracken_human_reads`,
`bracken_synthetic_reads` and their `pct_bracken_*` give the same two taxa from the Bracken
species table, as a share of all species-level reads.

| Parameter | Default | Meaning |
|---|---|---|
| `--consensus_min_reads` | `100000` | the MetaPhlAn gate above; a non-negative integer (`1e5` is rejected) |
| `--bracken_read_length` | `50` | Bracken `-r`; the Kraken database must hold a `database<N>mers.kmer_distrib` for it, checked at launch |
| `--bracken_threshold` | `2` | Bracken `-t`: the minimum reads a taxon needs to be re-estimated. A read threshold, not a thread count. A thin sample in which no taxon reaches it gets an empty no-call table instead of failing the run |

```bash
nextflow run main.nf --consensus_min_reads 100000 --bracken_read_length 50 --bracken_threshold 2
```

**Consensus** (`scripts/compute_consensus_taxa.py`) decides which genera are believable by
asking both profilers. It learns from the samples that have *both* a MetaPhlAn profile and a
Bracken table, which are exactly the samples that cleared `--consensus_min_reads`. Across those
samples, a genus is kept if it is present in at least 1% of them by Bracken **and** at least 1%
by MetaPhlAn. The consensus genus list is then applied to the Bracken tables of **all**
samples, including those below the gate: the output tables are Bracken counts restricted to
consensus genera (and, for species, to species of those genera). No taxon, including
*Cutibacterium*, is hard-coded out of the analysis.

**Pass-through.** With fewer than 2 such samples there is nothing to learn a consensus from,
so the Bracken tables pass through unfiltered, under the same file names, and the run continues.
This includes the case where no sample cleared the gate and MetaPhlAn never ran.
`CONSENSUS_TAXA/consensus_status.tsv` records which path ran, as counts only (no sample IDs):

| Key | Meaning |
|---|---|
| `mode` | `consensus` or `bracken_passthrough` |
| `reason` | why, when it passed through |
| `samples_in_bracken`, `samples_with_metaphlan`, `samples_used_for_consensus` | the sample counts at each step |
| `samples_above_min_reads` | samples whose Bracken **genus** total reaches the threshold. Reported for comparison only; the gate uses the KrakenUniq root count |
| `genus_rows_out`, `species_rows_out` | rows in the output tables |
| `species_without_genus_dropped` | non-viral Bracken species with no `g__` in their lineage (e.g. endosymbionts, `synthetic_construct`), left out because they cannot join a genus |

Check `mode` before interpreting a run's taxa: a pass-through table has not been filtered by
agreement between the two profilers. The cohort report shows the same status.

Consensus also writes `bracken_library_sizes.tsv`, each sample's total Bracken reads before any
filtering. Decontamination uses it to tell a sample that is all-zero after filtering from one
whose library was empty.

---
[Back to the README](../../README.md)
