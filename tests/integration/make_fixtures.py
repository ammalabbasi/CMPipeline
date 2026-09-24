#!/usr/bin/env python3
"""Build the synthetic fixtures for the executed integration test (audit C22).

Everything here is generated from a seeded random number generator: random "host" and
"microbe" genomes, reads simulated from them, a BAM and a CRAM built from those reads, paired
and single-end FASTQ, tiny databases, sample metadata and a decontamination input table. No
real sequence, sample, table or database is read or copied.

    python3 make_fixtures.py --outdir <fixture dir> [--seed 22]

Needs `samtools` and `minimap2` on PATH (run_integration.sh puts the pipeline's own envs there).
Writes <outdir>/expected.json with the exact counts check_outputs.py asserts on.
"""

import argparse
import gzip
import hashlib
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

READ_LEN = 100
FRAGMENT = 300
COMP = str.maketrans("ACGT", "TGCA")

# Synthetic taxonomy: (taxid, parent, KrakenUniq rank word, name). All names are invented.
TAXONOMY = [
    (1, 0, "no rank", "root"),
    (2, 1, "superkingdom", "Bacteria"),
    (1001, 2, "phylum", "Synthphylum"),
    (1002, 1001, "class", "Synthclass"),
    (1003, 1002, "order", "Synthorder"),
    (1004, 1003, "family", "Synthfamily"),
    (2001, 1004, "genus", "Alphagenus"),
    (2002, 1004, "genus", "Betagenus"),
    (2003, 1004, "genus", "Gammagenus"),
    (3001, 2001, "species", "Alphagenus prima"),
    (3002, 2001, "species", "Alphagenus secunda"),
    (3003, 2002, "species", "Betagenus tertia"),
    (3004, 2003, "species", "Gammagenus quarta"),
]
SPECIES = [3001, 3002, 3003, 3004]
GENUS_OF = {3001: 2001, 3002: 2001, 3003: 2002, 3004: 2003}

# Per-sample design. mic_pairs / mic_singles: microbial fragments by species taxid (a key
# "g2001" means reads the stub classifier assigns at GENUS level). host_*: host reads that are
# still present in the input and must be removed by the matching host index.
SAMPLES = {
    "itest_bam": dict(kind="bam", mic_pairs={3001: 10, 3002: 5, 3003: 4, "g2001": 2},
                      mic_singles={3001: 3, 3004: 2}, half_pair_taxon=3003,
                      host_a_unmapped_singles=2, host_b_unmapped_pairs=3, mapped_host_pairs=20),
    "itest_cram": dict(kind="cram", mic_pairs={3001: 6, 3003: 8, 3004: 4},
                       mic_singles={3002: 2}, half_pair_taxon=3001,
                       host_a_unmapped_singles=1, host_b_unmapped_pairs=2, mapped_host_pairs=15),
    "itest_fq1": dict(kind="fastq_paired", mic_pairs={3001: 12, 3002: 3, 3004: 6},
                      host_a_pairs=10, host_b_pairs=5),
    "itest_fq2": dict(kind="fastq_paired", mic_pairs={3001: 4, 3003: 10, 3004: 3},
                      host_a_pairs=8, host_b_pairs=4),
    # One microbial read: below Bracken's -t 2 at every rank, so Bracken says "no reads found"
    # and the module must write a no-call instead of failing the cohort (audit C02).
    "itest_thin": dict(kind="fastq_single", mic_singles={3004: 1}, host_a_singles=2),
}
# Decontam roles/batches for the full run (run A): each batch holds >= 1 positive and 1 control.
RUN_A_META = {
    "itest_bam": ("Tumor", "b1", 61, "F", 24.1),
    "itest_fq1": ("Blood", "b1", 55, "M", 27.9),
    "itest_thin": ("Tumor", "b1", 70, "M", 22.5),
    "itest_cram": ("Tumor", "b2", 48, "F", 30.2),
    "itest_fq2": ("Blood", "b2", 66, "F", 25.0),
}


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def run(cmd):
    subprocess.run(cmd, check=True)


class Simulator:
    def __init__(self, rng, genomes):
        self.rng = rng
        self.genomes = genomes      # name -> sequence
        self.serial = 0

    def _pos(self, genome, span):
        return self.rng.randrange(0, len(self.genomes[genome]) - span)

    def name(self, prefix):
        self.serial += 1
        return f"{prefix}_{self.serial}"

    def pair(self, genome, prefix):
        """Fragment -> (name, r1, r2, pos1, pos2): r1 forward, r2 reverse strand."""
        g = self.genomes[genome]
        p = self._pos(genome, FRAGMENT)
        r1 = g[p:p + READ_LEN]
        r2_fwd = g[p + FRAGMENT - READ_LEN:p + FRAGMENT]
        return self.name(prefix), r1, revcomp(r2_fwd), p, p + FRAGMENT - READ_LEN

    def single(self, genome, prefix):
        g = self.genomes[genome]
        p = self._pos(genome, READ_LEN)
        return self.name(prefix), g[p:p + READ_LEN], p


def mic_genome(key):
    # Genus-level reads are simulated from the genus's first species; only the name differs.
    if isinstance(key, str):
        genus = int(key[1:])
        return f"sp{min(s for s in SPECIES if GENUS_OF[s] == genus)}"
    return f"sp{key}"


def mic_prefix(key):
    return f"mic_{key}" if isinstance(key, str) else f"mic_s{key}"


QUAL = "I" * READ_LEN


def sam_line(name, flag, rname, pos, cigar, rnext, pnext, tlen, seq):
    mapq = 60 if not (flag & 4) else 0
    return "\t".join(map(str, [name, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, QUAL])) + "\n"


def build_alignment(sim, spec, host_seqs, prefix):
    """SAM text for one BAM/CRAM sample plus its expected counts."""
    lines = []
    exp = dict(primary=0, unmapped=0, host_a_unmapped=0, host_b_unmapped=0, microbial=0)
    # Mapped host pairs (flag 99/147), primary.
    for _ in range(spec["mapped_host_pairs"]):
        n, r1, r2, p1, p2 = sim.pair("hostA", "hstA")
        lines.append(sam_line(n, 99, "hostA", p1 + 1, "100M", "=", p2 + 1, FRAGMENT, r1))
        lines.append(sam_line(n, 147, "hostA", p2 + 1, "100M", "=", p1 + 1, -FRAGMENT, revcomp(r2)))
        exp["primary"] += 2
    # One secondary (256) and one supplementary (2048) record of a mapped host read: never
    # counted as library records and never extracted (-F 2304).
    n, r1, r2, p1, p2 = sim.pair("hostA", "hstA")
    lines.append(sam_line(n, 99, "hostA", p1 + 1, "100M", "=", p2 + 1, FRAGMENT, r1))
    lines.append(sam_line(n, 147, "hostA", p2 + 1, "100M", "=", p1 + 1, -FRAGMENT, revcomp(r2)))
    exp["primary"] += 2
    lines.append(sam_line(n, 99 + 256, "hostA", p1 + 1, "100M", "=", p2 + 1, FRAGMENT, r1))
    lines.append(sam_line(n, 147 + 2048, "hostA", p2 + 1, "60M40S", "=", p1 + 1, -FRAGMENT, revcomp(r2)))
    # Half-mapped pair: mate 1 on host (73), mate 2 unmapped microbial (133).
    n, r1, _, p1, _ = sim.pair("hostA", "hstA")
    p = p1 + 1
    _, m1, _ = sim.single(mic_genome(spec["half_pair_taxon"]), "x")
    half_name = n.replace("hstA", mic_prefix(spec["half_pair_taxon"]))
    lines.append(sam_line(half_name, 73, "hostA", p, "100M", "=", p, 0, r1))
    lines.append(sam_line(half_name, 133, "hostA", p, "*", "=", p, 0, m1))
    exp["primary"] += 2
    exp["unmapped"] += 1
    exp["microbial"] += 1
    # Unmapped microbial pairs (77/141).
    for key, count in spec["mic_pairs"].items():
        for _ in range(count):
            n, r1, r2, _, _ = sim.pair(mic_genome(key), mic_prefix(key))
            lines.append(sam_line(n, 77, "*", 0, "*", "*", 0, 0, r1))
            lines.append(sam_line(n, 141, "*", 0, "*", "*", 0, 0, r2))
            exp["primary"] += 2
            exp["unmapped"] += 2
            exp["microbial"] += 2
    # Flag-0-category unmapped reads (flag 4 only: neither READ1 nor READ2). The F01 defect
    # routed exactly these to stdout instead of the FASTQ.
    for key, count in spec["mic_singles"].items():
        for _ in range(count):
            n, s, _ = sim.single(mic_genome(key), mic_prefix(key))
            lines.append(sam_line(n, 4, "*", 0, "*", "*", 0, 0, s))
            exp["primary"] += 1
            exp["unmapped"] += 1
            exp["microbial"] += 1
    # Host reads the original aligner left unmapped: removed by the hg38 / t2t indexes.
    for _ in range(spec["host_a_unmapped_singles"]):
        n, s, _ = sim.single("hostA", "hstA")
        lines.append(sam_line(n, 4, "*", 0, "*", "*", 0, 0, s))
        exp["primary"] += 1
        exp["unmapped"] += 1
        exp["host_a_unmapped"] += 1
    for _ in range(spec["host_b_unmapped_pairs"]):
        n, r1, r2, _, _ = sim.pair("hostB", "hstB")
        lines.append(sam_line(n, 77, "*", 0, "*", "*", 0, 0, r1))
        lines.append(sam_line(n, 141, "*", 0, "*", "*", 0, 0, r2))
        exp["primary"] += 2
        exp["unmapped"] += 2
        exp["host_b_unmapped"] += 2
    return lines, exp


def write_fastq_gz(path, reads):
    with gzip.open(path, "wt") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{QUAL}\n")


def taxon_of_key(key):
    return int(key[1:]) if isinstance(key, str) else key


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=22)
    args = ap.parse_args()

    for tool in ("samtools", "minimap2"):
        if not shutil.which(tool):
            sys.exit(f"ERROR: {tool} not on PATH")

    out = Path(args.outdir).resolve()
    for sub in ("ref", "db/kraken", "db/metaphlan", "samples", "decontam", "batch"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    # ---- references: random bases only
    host = {"hostA": rand_seq(rng, 60000), "hostB": rand_seq(rng, 40000)}
    microbes = {f"sp{t}": rand_seq(rng, 20000) for t in SPECIES}
    write_fasta(out / "ref/host_a.fa", [("hostA", host["hostA"])])
    write_fasta(out / "ref/host_b.fa", [("hostB", host["hostB"])])
    write_fasta(out / "ref/host.fa", list(host.items()))       # CRAM reference: both contigs
    write_fasta(out / "ref/microbes.fa", list(microbes.items()))
    run(["samtools", "faidx", str(out / "ref/host.fa")])

    # ---- tiny host databases, built the way scripts/create_minimap2_indexes.sh builds them
    run(["minimap2", "-x", "sr", "-d", str(out / "db/hg38.mmi"), str(out / "ref/host_a.fa")])
    run(["minimap2", "-x", "sr", "-d", str(out / "db/t2t_phix.mmi"), str(out / "ref/host_b.fa")])

    # ---- Kraken/Bracken DB: taxonomy for the stub, a kmer distribution real Bracken can read
    kdb = out / "db/kraken"
    with open(kdb / "stub_taxonomy.tsv", "w") as fh:
        for row in TAXONOMY:
            fh.write("\t".join(map(str, row)) + "\n")
    (kdb / "taxDB").write_text("synthetic placeholder: only its existence is checked at launch\n")
    with open(kdb / "database50mers.kmer_distrib", "w") as fh:
        fh.write("mapped_taxid\tgenome_taxids:kmers_mapped:total_genome_kmers\n")
        for sp in SPECIES:
            fh.write(f"{sp}\t{sp}:1000:1000\n")
        for genus in sorted(set(GENUS_OF.values())):
            kids = [s for s in SPECIES if GENUS_OF[s] == genus]
            fh.write(f"{genus}\t" + " ".join(f"{s}:{1000 // len(kids)}:1000" for s in kids) + "\n")
    # MetaPhlAn is never run (no sample reaches --consensus_min_reads); launch checks the .pkl exists.
    (out / "db/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl").write_text("placeholder\n")

    # ---- samples
    genomes = dict(host)
    genomes.update(microbes)
    sim = Simulator(rng, genomes)
    expected = {"samples": {}}
    sheet = ["patient,bam,cram,cram_reference,fastq_1,fastq_2"]
    md5 = {k: hashlib.md5(v.encode()).hexdigest() for k, v in host.items()}
    header = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(
        f"@SQ\tSN:{k}\tLN:{len(v)}\tM5:{md5[k]}\n" for k, v in host.items())

    for sample, spec in SAMPLES.items():
        e = {"kind": spec["kind"]}
        taxa = {}
        if spec["kind"] in ("bam", "cram"):
            lines, c = build_alignment(sim, spec, host, sample)
            sam = out / f"samples/{sample}.sam"
            sam.write_text(header + "".join(lines))
            bam = out / f"samples/{sample}.bam"
            run(["samtools", "sort", "-o", str(bam), str(sam)])
            if spec["kind"] == "cram":
                cram = out / f"samples/{sample}.cram"
                run(["samtools", "view", "-C", "-T", str(out / "ref/host.fa"), "-o", str(cram), str(bam)])
                bam.unlink()
                sheet.append(f"{sample},,{cram},{out / 'ref/host.fa'},,")
            else:
                sheet.append(f"{sample},{bam},,,,")
            sam.unlink()
            e.update(library_primary_records=c["primary"], extracted_unmapped_records=c["unmapped"],
                     reads_raw=c["unmapped"],
                     reads_after_hg38=c["unmapped"] - c["host_a_unmapped"],
                     reads_after_t2t_phix=c["unmapped"] - c["host_a_unmapped"] - c["host_b_unmapped"],
                     reads_host_depleted=c["microbial"])
            for key, n in spec["mic_pairs"].items():
                taxa[taxon_of_key(key)] = taxa.get(taxon_of_key(key), 0) + 2 * n
            for key, n in spec["mic_singles"].items():
                taxa[taxon_of_key(key)] = taxa.get(taxon_of_key(key), 0) + n
            taxa[spec["half_pair_taxon"]] = taxa.get(spec["half_pair_taxon"], 0) + 1
        else:
            r1, r2 = [], []
            host_a = host_b = 0
            for key, n in spec.get("mic_pairs", {}).items():
                for _ in range(n):
                    name, a, b, _, _ = sim.pair(mic_genome(key), mic_prefix(key))
                    r1.append((name, a)); r2.append((name, b))
                taxa[taxon_of_key(key)] = taxa.get(taxon_of_key(key), 0) + 2 * n
            for _ in range(spec.get("host_a_pairs", 0)):
                name, a, b, _, _ = sim.pair("hostA", "hstA")
                r1.append((name, a)); r2.append((name, b)); host_a += 2
            for _ in range(spec.get("host_b_pairs", 0)):
                name, a, b, _, _ = sim.pair("hostB", "hstB")
                r1.append((name, a)); r2.append((name, b)); host_b += 2
            for key, n in spec.get("mic_singles", {}).items():
                for _ in range(n):
                    name, s, _ = sim.single(mic_genome(key), mic_prefix(key))
                    r1.append((name, s))
                taxa[taxon_of_key(key)] = taxa.get(taxon_of_key(key), 0) + n
            for _ in range(spec.get("host_a_singles", 0)):
                name, s, _ = sim.single("hostA", "hstA")
                r1.append((name, s)); host_a += 1
            f1 = out / f"samples/{sample}_R1.fastq.gz"
            write_fastq_gz(f1, r1)
            if r2:
                f2 = out / f"samples/{sample}_R2.fastq.gz"
                write_fastq_gz(f2, r2)
                sheet.append(f"{sample},,,,{f1},{f2}")
            else:
                sheet.append(f"{sample},,,,{f1},")
            total = len(r1) + len(r2)
            e.update(library_primary_records="NA", extracted_unmapped_records="NA", reads_raw=total,
                     reads_after_hg38=total - host_a, reads_after_t2t_phix=total - host_a - host_b,
                     reads_host_depleted=total - host_a - host_b)
        e["microbial_reads"] = sum(taxa.values())
        e["taxa"] = {str(k): v for k, v in sorted(taxa.items())}
        expected["samples"][sample] = e

    (out / "samplesheet.csv").write_text("\n".join(sheet) + "\n")
    with open(out / "metadata_runA.tsv", "w") as fh:
        fh.write("sampleid\tsample_type\tshipment_batch\tage_diag\tsex\tbmi\n")
        for s, (t, b, age, sex, bmi) in RUN_A_META.items():
            fh.write(f"{s}\t{t}\t{b}\t{age}\t{sex}\t{bmi}\n")

    # ---- run B: a standalone decontam input, 2 batches x 6 samples (4 positive, 2 control),
    # the shape consensus_taxa writes (clade_name + one integer column per sample).
    genera = [f"Synthgenus{c}" for c in "ABCDEFG"] + ["Contamgenus"]
    samples_b = [f"dc_{b}_{i}" for b in ("b1", "b2") for i in range(6)]
    roles_b = {s: ("Blood" if s.endswith(("_4", "_5")) else "Tumor") for s in samples_b}
    counts = {}
    for g in genera:
        for s in samples_b:
            if g == "Contamgenus":
                # Everywhere in controls, rarely and thinly in positives: the prevalence signal.
                v = rng.randint(40, 80) if roles_b[s] == "Blood" else (rng.randint(5, 9) if rng.random() < 0.25 else 0)
            else:
                v = rng.randint(20, 400) if roles_b[s] == "Tumor" or rng.random() < 0.3 else 0
            counts[(g, s)] = v
    lineage = "d__Bacteria|p__Synthphylum|c__Synthclass|o__Synthorder|f__Synthfamily|g__"
    with open(out / "decontam/genus_table.tsv", "w") as fh:
        fh.write("clade_name\t" + "\t".join(samples_b) + "\n")
        for g in genera:
            fh.write(lineage + g + "\t" + "\t".join(str(counts[(g, s)]) for s in samples_b) + "\n")
    with open(out / "decontam/library_sizes.tsv", "w") as fh:
        fh.write("sample\tbracken_genus_total\tbracken_species_total\n")
        for s in samples_b:
            tot = sum(counts[(g, s)] for g in genera) + 5000
            fh.write(f"{s}\t{tot}\t{tot}\n")
    with open(out / "decontam/metadata.tsv", "w") as fh:
        fh.write("sampleid\tsample_type\tshipment_batch\tage_diag\tsex\tbmi\n")
        for i, s in enumerate(samples_b):
            fh.write(f"{s}\t{roles_b[s]}\t{s.split('_')[1]}\t{40 + 3 * i}\t{'FM'[i % 2]}\t{21.0 + 0.7 * i:.1f}\n")
    expected["decontam_run_b"] = {"samples": len(samples_b), "batches": 2, "taxa": len(genera)}

    # ---- run C (audit C22): enough samples for ConQuR (>= 20 after decontam) -- 2 batches x 14
    # (12 positive, 2 control) -- with a real batch shift (batch 2 counts x3 on half the genera)
    # and covariates that vary within each batch, so batch correction must actually CORRECT.
    samples_c = [f"bc_{b}_{i}" for b in ("b1", "b2") for i in range(14)]
    roles_c = {s: ("Blood" if s.endswith(("_12", "_13")) else "Tumor") for s in samples_c}
    counts_c = {}
    for j, g in enumerate(genera[:-1]):
        for s in samples_c:
            v = rng.randint(50, 500) if roles_c[s] == "Tumor" or rng.random() < 0.3 else 0
            if s.startswith("bc_b2") and j % 2 == 0:
                v *= 3
            counts_c[(g, s)] = v
    with open(out / "batch/genus_table.tsv", "w") as fh:
        fh.write("clade_name\t" + "\t".join(samples_c) + "\n")
        for g in genera[:-1]:
            fh.write(lineage + g + "\t" + "\t".join(str(counts_c[(g, s)]) for s in samples_c) + "\n")
    with open(out / "batch/library_sizes.tsv", "w") as fh:
        fh.write("sample\tbracken_genus_total\tbracken_species_total\n")
        for s in samples_c:
            tot = sum(counts_c[(g, s)] for g in genera[:-1]) + 5000
            fh.write(f"{s}\t{tot}\t{tot}\n")
    for name, cols in (("metadata.tsv", ("age_diag", "sex", "bmi")), ("metadata_nocov.tsv", ("age_diag", "sex"))):
        with open(out / "batch" / name, "w") as fh:
            fh.write("sampleid\tsample_type\tshipment_batch\t" + "\t".join(cols) + "\n")
            for i, s in enumerate(samples_c):
                vals = {"age_diag": str(35 + 2 * i), "sex": "FM"[i % 2], "bmi": f"{20.0 + 0.5 * i:.1f}"}
                fh.write(f"{s}\t{roles_c[s]}\t{s.split('_')[1]}\t" + "\t".join(vals[c] for c in cols) + "\n")
    expected["batch_run_c"] = {"samples": len(samples_c), "batches": 2}

    (out / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")
    print(f"fixtures: {len(SAMPLES)} pipeline samples, {len(samples_b)} decontam samples -> {out}")


if __name__ == "__main__":
    main()
