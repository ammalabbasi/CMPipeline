# Classification and consensus

How KrakenUniq/Bracken, the MetaPhlAn gate and consensus behave, plus adapter trimming and empty samples.

## Adapter trimming

`--adapter_trim` controls whether fastp scans reads against `--adapters`
(`ref/known_adapters.fna`, 234 sequences):

- `auto` (default) scans raw FASTQ input, and skips the scan for reads extracted
  from a BAM or CRAM, which were adapter-trimmed before they were aligned.
- `always` scans every sample; `never` scans none.

The scan is expensive and, on aligned input, produces nothing: measured on one
aligned CRAM, `filterReads` took 2 h 50 m with `--adapter_fasta` and 3 m 15 s
without it, with fastp reporting no adapter detected and 0 reads trimmed.
Skipping the FASTA scan does not disable adapter trimming; fastp still
auto-detects adapters.

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
`bracken_microbial_reads` and `cleared_consensus_gate` columns show, per sample, the count and
whether it cleared.

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

Check `mode` before interpreting a run's taxa: a pass-through table has not been filtered by
agreement between the two profilers. The cohort report shows the same status.

Consensus also writes `bracken_library_sizes.tsv`, each sample's total Bracken reads before any
filtering. Decontamination uses it to tell a sample that is all-zero after filtering from one
whose library was empty.

---
[Back to the README](../../README.md)
