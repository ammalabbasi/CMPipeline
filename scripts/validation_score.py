#!/usr/bin/env python3
"""Alignment validation: score one sample's alignments against the competitive reference.

Reads minimap2 SAM on stdin (or --sam), aligned with
`minimap2 -ax sr -p 0.99 -N 100 --secondary=yes` against the validation index, and writes
<out_dir>/<sample>.validation_stats.tsv: one row per taxon that sample's Bracken table reports.

One sample (--sample S): read names are used as they are.
Several samples in ONE stream (--samples S1,S2,...): read names must be `<sample>::<read>`, with
each sample's reads contiguous (samples are concatenated in order). This lets one minimap2 run
serve a whole batch, so the large index is loaded once per batch instead of once per sample.
Every listed sample gets a stats file, even one with no reads at all.

Per read (a read = one primary record plus its secondaries; mates of a pair are separate reads,
which is how host depletion hands them over, as one interleaved FASTQ):
  best hits = alignments with the read's best AS; a taxon is "supported" by the read if any best
  hit is on it; the read is "unique" to a taxon if every best hit is on that taxon. Taxa in one
  --groups group (default: E. coli + Shigella) count as one taxon throughout.

Per taxon, measured on its best-supported genome (the accession with most supporting reads):
  breadth_ratio    covered bases / expected under uniform sampling, 1 - exp(-N*L/G)
  bins_ratio       10 kb windows hit / expected, counting each FRAGMENT (read name) once, since
                   mates ~500 bp apart would otherwise halve the apparent spread
  median_identity  1 - NM / aligned columns, over all supporting alignments
  unique_fraction  unique / (unique + shared) reads

Columns: supporting_reads, unique_reads, multi_reads, unique_fraction and median_identity POOL every
genome of the taxon; breadth_ratio, bins_ratio, reads_on_best and genome_len use only the best
genome. minimap2_count equals unique_reads (reported, never used). (audit B08)

Status, tested in this order (thresholds are arguments; defaults from the 2026-09-23 benchmark):
  unsupported       Bracken >= --min_support but supporting reads < --min_confirmed (10): the reads are
                    not this species'; zeroed like `failed` (decision 29, as PRISM)
  insufficient      supporting reads over all genomes < --min_support (Bracken too low to test, or
                    --min_confirmed..--min_support-1 supporting reads)
  failed            unique_fraction < --min_unique_frac (a shadow of another species; audit B03)
  insufficient      reads on the best genome < --min_support (spread thin over many genomes)
  validated         bins_ratio >= --min_bins_ratio and identity >= --min_identity
  failed            otherwise; fail_reason lists every test that failed
  not_in_reference  the taxon has no genome in the validation reference, or is Homo_sapiens

With --export_reads, reads supporting a VALIDATED taxon are written (original orientation) to a
gzipped FASTA, headers  >read_name taxon=<label>. These are sample-derived reads: protected data.
"""
import argparse, collections, gzip, math, os, re, shutil, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_common import (HUMAN_LABEL, label_of, read_groups, read_species_table,
                               species_from_clade)

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def ref_len(cigar):
    return sum(int(n) for n, op in CIGAR.findall(cigar) if op in "MDN=X")


def aln_columns(cigar):
    return sum(int(n) for n, op in CIGAR.findall(cigar) if op in "MID=X")


def merged_length(intervals):
    total, cs, ce = 0, None, None
    for s, e in sorted(intervals):
        if ce is None or s > ce:
            if ce is not None:
                total += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return total + (ce - cs if ce is not None else 0)


def bins_expected(contig_lens, G, n, bin_size):
    e = 0.0
    for clen in contig_lens:
        full, rest = divmod(clen, bin_size)
        for w in [bin_size] * full + ([rest] if rest else []):
            e += 1 - (1 - w / G) ** n
    return e


def median(v):
    v = sorted(v)
    k = len(v)
    if not k:
        return float("nan")
    return v[k // 2] if k % 2 else (v[k // 2 - 1] + v[k // 2]) / 2


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sam", default="-")
    ap.add_argument("--sample", default=None, help="single sample; read names unprefixed")
    ap.add_argument("--samples", default=None, help="comma-separated; read names <sample>::<read>")
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--contigs", required=True, help="contig_labels.tsv from the reference build")
    ap.add_argument("--species_table", required=True, help="the species table going downstream")
    ap.add_argument("--groups", default=None)
    ap.add_argument("--min_support", type=int, default=25)
    ap.add_argument("--min_confirmed", type=int, default=10,
                    help="with >= min_support Bracken reads, fewer supporting reads than this is 'unsupported' (zeroed)")
    ap.add_argument("--min_bins_ratio", type=float, default=0.5)
    ap.add_argument("--min_identity", type=float, default=0.98)
    ap.add_argument("--min_unique_frac", type=float, default=0.02)
    ap.add_argument("--bin_size", type=int, default=10_000)
    ap.add_argument("--export_reads", action="store_true",
                    help="write <sample>.validated_reads.fasta.gz (reads supporting validated taxa)")
    a = ap.parse_args(argv)
    if bool(a.sample) == bool(a.samples):
        sys.exit("ERROR: give exactly one of --sample or --samples")

    groups = read_groups(a.groups)
    def taxon(label):
        return groups.get(label, label)

    genome_len = collections.Counter()
    contig_lens = collections.defaultdict(list)
    ref_taxa = set()
    with open(a.contigs) as fh:
        fh.readline()
        for line in fh:
            contig, label, acc, n = line.rstrip("\n").split("\t")[:4]
            genome_len[acc] += int(n)
            contig_lens[acc].append(int(n))
            ref_taxa.add(taxon(label))

    table_samples, table = read_species_table(a.species_table)
    wanted = [a.sample] if a.sample else [x for x in a.samples.split(",") if x]
    absent = [x for x in wanted if x not in table_samples]
    if absent:
        sys.exit(f"ERROR: sample(s) not columns of {a.species_table}: {', '.join(absent)}")

    scorers = {x: SampleScorer(x, a, taxon, table, genome_len, contig_lens, ref_taxa) for x in wanted}
    done = set()
    current = None
    sam = sys.stdin if a.sam == "-" else open(a.sam)
    cur, alns = None, []

    def flush():
        nonlocal current
        if cur is None:
            return
        qname = cur[0]
        if a.samples:
            smp, sep, _ = qname.partition("::")
            if not sep or smp not in scorers:
                sys.exit(f"ERROR: read '{qname}' has no known <sample>:: prefix")
        else:
            smp = a.sample
        if smp != current:
            if smp in done:
                sys.exit(f"ERROR: reads of sample '{smp}' are not contiguous in the SAM stream")
            if current is not None:
                scorers[current].finish(); done.add(current)
            current = smp
        scorers[smp].add_read(qname, cur[1], cur[2], alns)

    for line in sam:
        if line.startswith("@"):
            continue
        f = line.split("\t", 11)
        flag = int(f[1])
        if flag & 2048:                    # supplementary: part of the same read, ignore
            continue
        if not flag & 256:                 # primary or unmapped record starts a new read
            flush()
            cur, alns = (f[0], f[9], bool(flag & 16)), []
            if flag & 4:
                continue
        m = re.search(r"\tAS:i:(-?\d+)", line)
        nm = re.search(r"\tNM:i:(\d+)", line)
        cols = aln_columns(f[5])
        ident = 1 - int(nm.group(1)) / cols if nm and cols else float("nan")
        alns.append((int(m.group(1)) if m else 0, f[2], int(f[3]) - 1, ref_len(f[5]),
                     not flag & 256, ident))
    flush()
    for x in wanted:                       # includes samples with no reads at all
        if x not in done:
            scorers[x].finish(); done.add(x)


# minimap2_count = reads minimap2 assigns to this taxon ALONE (= unique_reads): the alignment-based
# count the DESIGN promised, reported for sensitivity analysis and never used downstream (B08).
# supporting/unique/multi pool every genome of the taxon; breadth, bins and reads_on_best use the
# best-supported genome only -- as the column names say.
COLS = ["sample", "taxon", "members", "bracken_reads", "minimap2_count", "supporting_reads", "unique_reads",
        "multi_reads", "unique_fraction", "best_genome", "reads_on_best", "genome_len",
        "observed_breadth", "expected_breadth", "breadth_ratio", "bins_hit", "bins_expected",
        "bins_ratio", "median_identity", "status", "fail_reason"]


class SampleScorer:
    """Accumulates one sample's reads; finish() writes its stats (and exported reads)."""

    def __init__(self, sample, a, taxon, table, genome_len, contig_lens, ref_taxa):
        self.sample, self.a, self.taxon = sample, a, taxon
        self.genome_len, self.contig_lens, self.ref_taxa = genome_len, contig_lens, ref_taxa
        self.bracken = collections.Counter(); self.members = collections.defaultdict(set)
        for clade, counts in table.items():
            lab = label_of(species_from_clade(clade))
            self.bracken[taxon(lab)] += counts[sample]
            self.members[taxon(lab)].add(lab)
        self.uniq = collections.Counter(); self.multi = collections.Counter()
        # (taxon, accession) -> list of (contig, start, end, identity, fragment_name, read_len)
        self.support = collections.defaultdict(list)
        self.tmpdir = tempfile.mkdtemp(prefix="valreads_") if a.export_reads else None
        self.handles = {}

    def add_read(self, qname, seq, rev, alns):
        if not alns:
            return
        taxon = self.taxon
        best = max(x[0] for x in alns)
        top = [x for x in alns if x[0] == best]
        taxa = {taxon(x[1].split("|", 1)[0]) for x in top}
        for t in taxa:
            hit = next((x for x in top if taxon(x[1].split("|", 1)[0]) == t and x[4]),
                       next(x for x in top if taxon(x[1].split("|", 1)[0]) == t))
            acc = hit[1].split("|")[1]
            self.support[(t, acc)].append((hit[1], hit[2], hit[2] + hit[3], hit[5], qname, hit[3]))
            # Export only reads that are NOT tied with human (audit B10): a read whose best hits
            # include human is not evidence for a microbe, whatever else it also matches.
            if self.tmpdir and t != HUMAN_LABEL and HUMAN_LABEL not in taxa and seq and seq != "*":
                if t not in self.handles:   # temp files named by index: labels may hold '/'
                    self.handles[t] = open(os.path.join(self.tmpdir, f"taxon_{len(self.handles)}.fa"), "w")
                s = seq.translate(COMP)[::-1] if rev else seq
                self.handles[t].write(f">{qname.partition('::')[2] or qname} taxon={t}\n{s}\n")
        if len(taxa) == 1:
            self.uniq[next(iter(taxa))] += 1
        else:
            for t in taxa:
                self.multi[t] += 1

    def finish(self):
        a, uniq, multi, support = self.a, self.uniq, self.multi, self.support
        for h in self.handles.values():
            h.close()
        validated = []
        with open(os.path.join(a.out_dir, f"{self.sample}.validation_stats.tsv"), "w") as out:
            out.write("\t".join(COLS) + "\n")
            for t in sorted(self.bracken, key=lambda k: -self.bracken[k]):
                if self.bracken[t] <= 0:
                    continue
                row = dict(sample=self.sample, taxon=t, members=",".join(sorted(self.members[t])),
                           bracken_reads=f"{self.bracken[t]:g}", minimap2_count=uniq[t],
                           unique_reads=uniq[t], multi_reads=multi[t])
                sup = [(len(v), acc) for (tt, acc), v in support.items() if tt == t]
                row["supporting_reads"] = sum(n for n, _ in sup)
                tot = uniq[t] + multi[t]
                row["unique_fraction"] = f"{uniq[t] / tot:.4f}" if tot else "NA"
                fail = []
                # Human is in the index only as the competitor that soaks up host reads; a Homo_sapiens
                # row in a pass-through / --skip_consensus table is never scored or masked (audit B09).
                # Bracken reported enough reads to test (>= min_support), yet alignment backs almost
                # none of them: the species is not supported, and its count is zeroed like `failed`
                # (decision 29; PRISM excludes taxa with < 10 uniquely confirmed reads).
                tested = self.bracken[t] >= a.min_support
                if t not in self.ref_taxa or t == HUMAN_LABEL:
                    status = "not_in_reference"
                elif not sup:
                    status = "unsupported" if tested else "insufficient"
                    if tested:
                        fail.append(f"supporting_reads<{a.min_confirmed}")
                else:
                    n, acc = max(sup)
                    hits = support[(t, acc)]
                    by_contig = collections.defaultdict(list); frag_bin = {}
                    for c, st, en, ident, q, rl in hits:
                        by_contig[c].append((st, en))
                        frag_bin.setdefault(q.rsplit("/", 1)[0], (c, st // a.bin_size))
                    G = self.genome_len[acc]
                    L = sum(h[5] for h in hits) / len(hits)
                    obs = sum(merged_length(v) for v in by_contig.values()) / G
                    exp = 1 - math.exp(-n * L / G)
                    be = bins_expected(self.contig_lens[acc], G, len(frag_bin), a.bin_size)
                    bins_ratio = len(set(frag_bin.values())) / be if be else float("nan")
                    med = median([h[3] for (tt, _), v in support.items() if tt == t for h in v])
                    ufrac = uniq[t] / tot if tot else 0.0
                    row.update(best_genome=acc, reads_on_best=n, genome_len=G,
                               observed_breadth=f"{obs:.6f}", expected_breadth=f"{exp:.6f}",
                               breadth_ratio=f"{obs / exp:.4f}" if exp else "NA",
                               bins_hit=len(set(frag_bin.values())), bins_expected=f"{be:.2f}",
                               bins_ratio=f"{bins_ratio:.4f}", median_identity=f"{med:.4f}")
                    total = row["supporting_reads"]
                    # Order matters (audit B03, 2026-09-23). A shadow -- every read shared with
                    # another taxon -- is detectable from all its supporting reads, however
                    # thinly they spread over its many strain genomes. Spread and identity need
                    # depth on ONE genome, gated exactly as in the benchmark (reads on best).
                    if total < a.min_support:
                        if tested and total < a.min_confirmed:
                            status = "unsupported"
                            fail.append(f"supporting_reads<{a.min_confirmed}")
                        else:
                            status = "insufficient"
                    elif not ufrac >= a.min_unique_frac:
                        fail.append(f"unique_fraction<{a.min_unique_frac:g}")
                        status = "failed"
                    elif n < a.min_support:
                        status = "insufficient"
                    else:
                        if not bins_ratio >= a.min_bins_ratio:
                            fail.append(f"bins_ratio<{a.min_bins_ratio:g}")
                        if not med >= a.min_identity:
                            fail.append(f"identity<{a.min_identity:g}")
                        status = "failed" if fail else "validated"
                if status == "validated":
                    validated.append(t)
                row["status"] = status
                row["fail_reason"] = ";".join(fail)
                out.write("\t".join(str(row.get(c, "NA")) for c in COLS) + "\n")
        if self.tmpdir:
            with gzip.open(os.path.join(a.out_dir, f"{self.sample}.validated_reads.fasta.gz"), "wt") as out:
                for t in validated:
                    if t in self.handles:
                        with open(self.handles[t].name) as fh:
                            shutil.copyfileobj(fh, out)
            shutil.rmtree(self.tmpdir)
        self.support.clear()                  # free memory before the next sample


if __name__ == "__main__":
    main()
