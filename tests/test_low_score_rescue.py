#!/usr/bin/env python3
"""Low-score mapped rescue (--extract_low_score_mapped auto|on|off; Ludmil 1 Oct, L8).

GDC aligned TCGA with `bwa mem -T 0`. bwa reports a read unmapped only when its best local score
is below -T, so under -T 0 microbial reads with a chance short human hit come out MAPPED (AS < 30)
and FLAG 4 extraction never sees them. The rescue adds PRIMARY mapped records with AS < 30 -- the
set a default bwa mem (-T 30) would have left unmapped -- and `auto` turns it on only for bwa mem /
bwa-mem2 mem with -T < 30. Ported from pksProfiler feature/e1-extraction (6138519, 43eac77).

The extraction tests render Modules/extract_reads.nf's whole script block exactly as Nextflow
would (Groovy escapes undone, the interpolated values substituted) and run it with real samtools
on a BAM and a CRAM built here from hand-written SAM records. Every record is synthetic: invented
names, random bases, one invented contig. Set CMP_TEST_SAMTOOLS_BIN to the pipeline's samtools env
bin/ (samtools, htsfile, bgzip, python) or have samtools >= 1.16 on PATH; otherwise they skip.
"""

import csv
import gzip
import hashlib
import importlib.util
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = (ROOT / "Modules" / "extract_reads.nf").read_text()
SURVIVAL_MODULE = (ROOT / "Modules" / "low_as_survival.nf").read_text()
MAIN = (ROOT / "main.nf").read_text()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


LSR = load("low_score_rescue", ROOT / "scripts" / "low_score_rescue.py")
QC = load("build_qc_summary_lsr", ROOT / "scripts" / "build_qc_summary.py")
FP = load("finalize_provenance_lsr", ROOT / "scripts" / "finalize_provenance.py")


def find_samtools_bin():
    candidates = [os.environ.get("CMP_TEST_SAMTOOLS_BIN"), str(Path(sys.executable).parent)]
    which = shutil.which("samtools")
    if which:
        candidates.append(str(Path(which).parent))
    for directory in candidates:
        if not directory:
            continue
        d = Path(directory)
        if not all((d / tool).is_file() for tool in ("samtools", "htsfile", "bgzip")):
            continue
        out = subprocess.run([str(d / "samtools"), "--version"], capture_output=True, text=True).stdout
        m = re.match(r"samtools (\d+)\.(\d+)", out)
        if m and (int(m.group(1)), int(m.group(2))) >= (1, 16):
            return d
    return None


SAMTOOLS_BIN = find_samtools_bin()

# ------------------------------------------------------------------ synthetic alignments

COMP = str.maketrans("ACGT", "TGCA")
RNG = random.Random(30)
CONTIG = "".join(RNG.choice("ACGT") for _ in range(5000))
QUAL = "I" * 59 + "#"   # asymmetric, so a reversed QUAL is visible


def revcomp(s):
    return s.translate(COMP)[::-1]


def seq(i, n=60):
    return "".join(random.Random(1000 + i).choice("ACGT") for _ in range(n))


# (name, flag, mapped?, AS or None, original read sequence, expected class)
# flag bits: 1 paired, 4 unmapped, 16 reverse, 32 mate reverse, 64 first, 128 second,
# 256 secondary, 2048 supplementary
RECORDS = [
    ("unm_pair", 77, False, None, seq(1), "flag4"),            # both unmapped, /1
    ("unm_pair", 141, False, None, seq(2), "flag4"),           # /2
    ("as29", 65 + 32, True, 29, seq(3), "low"),                # AS 29: rescued, /1
    ("as29", 145, True, 0, seq(4), "low"),                     # AS 0, reverse strand: rescued, /2
    ("as30", 99, True, 30, seq(5), "excluded"),                # AS 30: NOT below 30
    ("high", 65, True, 58, seq(6), "excluded"),                # AS 58
    ("no_as", 65, True, None, seq(7), "no_as"),                # mapped, no AS: never, counted
    ("secondary", 65 + 256, True, 5, seq(8), "excluded"),      # secondary, AS 5: never
    ("supplementary", 65 + 2048, True, 5, seq(9), "excluded"), # supplementary, AS 5: never
    ("unm_mate_mapped", 133 + 32, False, None, seq(10), "flag4"),  # FLAG 4 with a mapped mate
    ("unpaired_low", 16, True, 12, seq(11), "low"),            # unpaired, reverse, AS 12: no suffix
    ("unpaired_unm", 4, False, None, seq(12), "flag4"),        # unpaired FLAG 4: category 0
    ("unm_with_as", 4, False, 0, seq(13), "flag4"),            # FLAG 4 carrying AS: still FLAG 4
]
PRIMARY = sum(1 for _, flag, *_ in RECORDS if not flag & 2304)

PG_GDC = ("@PG\tID:bwa\tPN:bwa\tVN:0.7.15-r1140\tCL:bwa mem -t 36 -T 0 -R "
          "@RG\\tID:SECRET_RG\\tSM:SECRET_SAMPLE /ref/GRCh38.d1.vd1.fa SECRET_1.fq SECRET_2.fq")
PG_BWA_DEFAULT = "@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188\tCL:bwa mem -t 8 -K 100000000 ref.fa SECRET_a.fq SECRET_b.fq"
PG_BWA_T30 = "@PG\tID:bwa_mem\tPN:bwa\tVN:0.7.17\tCL:/opt/bin/bwa mem -Y -t 8 -T30 ref.fa in.fq"
PG_BWA_T19 = "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem -T 19 ref.fa in.fq"
PG_BWA_MEM2 = "@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:2.2.1\tCL:bwa-mem2 mem -t 8 -T 0 ref.fa a.fq b.fq"
PG_SENTIEON = "@PG\tID:x\tPN:sentieon\tVN:202112\tCL:sentieon bwa mem -T 0 r.fa a.fq"
PG_STAR = "@PG\tID:STAR\tPN:STAR\tVN:2.7.5c\tCL:STAR --runThreadN 8 --genomeDir /ref --readFilesIn SECRET.fq"
PG_MM2 = "@PG\tID:minimap2\tPN:minimap2\tVN:2.28-r1209\tCL:minimap2 -ax sr ref.mmi a.fq"
PG_ALN = "@PG\tID:bwa\tPN:bwa\tVN:0.7.12\tCL:bwa sampe ref.fa a.sai b.sai a.fq b.fq"
PG_NO_CL = "@PG\tID:bwa\tPN:bwa\tVN:0.7.17"
PG_SAMTOOLS = "@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.3\tCL:samtools sort -o SECRET.bam"
PG_UNKNOWN = "@PG\tID:mytool\tPN:mytool\tVN:1"


def header(*pg):
    md5 = hashlib.md5(CONTIG.encode()).hexdigest()
    return "\n".join(["@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:synthetic_chr1\tLN:{len(CONTIG)}\tM5:{md5}"]
                     + list(pg)) + "\n"


def sam_text(pg=(), records=RECORDS):
    lines = [header(*pg).rstrip("\n")]
    for name, flag, mapped, score, original, _ in records:
        stored = revcomp(original) if flag & 16 else original
        qual = QUAL[::-1] if flag & 16 else QUAL
        rname, pos, cigar, mapq = ("synthetic_chr1", 100, "20M40S", 0) if mapped else ("*", 0, "*", 0)
        if not mapped and flag & 1 and not flag & 8:
            rname, pos = "synthetic_chr1", 100   # placed at its mapped mate, as bwa does
        fields = [name, str(flag), rname, str(pos), str(mapq), cigar, "*", "0", "0", stored, qual]
        if score is not None:
            fields += ["NM:i:0", f"AS:i:{score}", "XS:i:0"]
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


def suffix(flag):
    return "/1" if flag & 64 and not flag & 128 else "/2" if flag & 128 and not flag & 64 else ""


def expected_reads(classes):
    """What samtools fastq -N writes: /1 /2 from the mate bits, reverse strand undone."""
    return sorted((name + suffix(flag), original, QUAL)
                  for name, flag, _m, _s, original, cls in RECORDS if cls in classes)


# ------------------------------------------------------------------ render the module script


def groovy_unescape(text):
    table = {"\\": "\\", "$": "$", "t": "\t", "n": "\n", '"': '"', "'": "'"}
    return re.sub(r"\\(.)", lambda m: table.get(m.group(1), m.group(0)), text)


def render_script(mode, alignment, reference="", input_type="bam", cpus=1):
    body = MODULE[MODULE.index('    """\n    set -euo pipefail') + len('    """\n'):]
    body = body[:body.index('    """\n}')]
    values = {
        "${sampleID}": "S", "${inputType}": input_type, "${alignment}": alignment,
        "${reference ?: ''}": reference, "${task.cpus}": str(cpus),
        "${params.scripts}": str(ROOT / "scripts"), "${params.extract_low_score_mapped}": mode,
    }
    for groovy, value in values.items():
        body = body.replace(groovy, value)
    left = re.findall(r"(?<!\\)\$\{[^}]*\}", body)
    assert not left, f"unsubstituted Groovy interpolation: {left}"
    return groovy_unescape(body)


def read_fastq_gz(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return sorted((lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4))


def read_tsv_row(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1, rows
    return rows[0]


@unittest.skipUnless(SAMTOOLS_BIN, "samtools >= 1.16 + htsfile + bgzip not found (set CMP_TEST_SAMTOOLS_BIN)")
class ExtractionOnSyntheticAlignments(unittest.TestCase):
    """The rendered extractReads script, run for real on synthetic BAM / CRAM."""

    def build(self, tmp, pg, kind="bam"):
        (tmp / "in.sam").write_text(sam_text(pg))
        tool = str(SAMTOOLS_BIN / "samtools")
        if kind == "bam":
            subprocess.run([tool, "view", "-b", "-o", str(tmp / "in.bam"), str(tmp / "in.sam")], check=True)
            return tmp / "in.bam", ""
        (tmp / "ref.fa").write_text(">synthetic_chr1\n" + CONTIG + "\n")
        subprocess.run([tool, "faidx", str(tmp / "ref.fa")], check=True)
        subprocess.run([tool, "view", "-C", "-T", str(tmp / "ref.fa"), "-o", str(tmp / "in.cram"),
                        str(tmp / "in.sam")], check=True)
        return tmp / "in.cram", str(tmp / "ref.fa")

    def run_module(self, mode, pg, kind="bam"):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        alignment, reference = self.build(tmp, pg, kind)
        work = tmp / "task"
        work.mkdir()
        (work / ".command.sh").write_text(render_script(mode, str(alignment), reference, kind))
        env = {"PATH": f"{SAMTOOLS_BIN}:/usr/bin:/bin", "HOME": str(tmp), "LC_ALL": "C"}
        result = subprocess.run(["bash", ".command.sh"], cwd=work, env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr[-2000:])
        self.assertEqual(result.stdout, "", "nothing may be written to the task's stdout (.command.out)")
        out = {
            "reads": read_fastq_gz(work / "S.UNMAPPED.fastq.gz"),
            "counts": read_tsv_row(work / "S.extract_counts.tsv"),
            "aligner": read_tsv_row(work / "S.upstream_aligner.tsv"),
            "stderr": result.stderr,
            "ids": None,
            "leftovers": sorted(p.name for p in work.iterdir()),
        }
        ids = work / "S.low_as_read_ids.txt.gz"
        if ids.exists():
            with gzip.open(ids, "rt") as fh:
                out["ids"] = sorted(fh.read().split())
        return out

    # --- the rule

    def test_auto_on_bwa_mem_T0_takes_flag4_and_primary_as_below_30(self):
        out = self.run_module("auto", (PG_GDC, PG_SAMTOOLS))
        self.assertEqual(out["reads"], expected_reads(("flag4", "low")))
        self.assertEqual(out["counts"], {
            "sample": "S", "library_primary_records": str(PRIMARY),
            "extracted_unmapped_records": "5", "low_as_reads_extracted": "3", "low_as_no_as_records": "1"})
        self.assertEqual(out["aligner"]["low_as_mode"], "on")
        self.assertEqual(out["aligner"]["low_as_mode_reason"], "auto_bwa_mem_T_below_30")
        self.assertIn("1 mapped primary records had no AS tag and were not added", out["stderr"])

    def test_as29_rescued_as30_not(self):
        names = {n.split("/")[0] for n, _, _ in self.run_module("on", (PG_GDC,))["reads"]}
        self.assertIn("as29", names)
        self.assertNotIn("as30", names)
        self.assertNotIn("high", names)

    def test_secondary_supplementary_and_missing_as_never_rescued(self):
        out = self.run_module("on", (PG_GDC,))
        names = {n.split("/")[0] for n, _, _ in out["reads"]}
        for never in ("secondary", "supplementary", "no_as"):
            self.assertNotIn(never, names)
        self.assertEqual(out["counts"]["low_as_no_as_records"], "1")

    def test_flag4_always_extracted(self):
        for mode, pg in (("auto", (PG_GDC,)), ("auto", (PG_STAR,)), ("off", (PG_GDC,)), ("on", ())):
            with self.subTest(mode=mode, pg=pg):
                reads = self.run_module(mode, pg)["reads"]
                self.assertTrue(set(expected_reads(("flag4",))) <= set(reads))

    def test_names_and_orientation_as_samtools_fastq_minus_n(self):
        reads = {n: (s, q) for n, s, q in self.run_module("on", (PG_GDC,))["reads"]}
        self.assertEqual(reads["as29/1"], (seq(3), QUAL))
        self.assertEqual(reads["as29/2"], (seq(4), QUAL))          # reverse strand undone
        self.assertEqual(reads["unpaired_low"], (seq(11), QUAL))   # no mate bit, no suffix

    def test_id_list_is_exactly_the_rescued_fastq_names(self):
        out = self.run_module("on", (PG_GDC,))
        self.assertEqual(out["ids"], sorted(n for n, _, _ in expected_reads(("low",))))
        self.assertTrue(set(out["ids"]) <= {n for n, _, _ in out["reads"]})

    def test_off_extracts_exactly_flag4_and_writes_no_id_list(self):
        for mode, pg in (("auto", (PG_STAR,)), ("auto", (PG_BWA_DEFAULT,)), ("auto", ()), ("off", (PG_GDC,))):
            with self.subTest(mode=mode, pg=pg):
                out = self.run_module(mode, pg)
                self.assertEqual(out["reads"], expected_reads(("flag4",)))
                self.assertEqual(out["counts"]["extracted_unmapped_records"], "5")
                self.assertEqual(out["counts"]["low_as_reads_extracted"], "0")
                self.assertEqual(out["counts"]["low_as_no_as_records"], "NA")
                self.assertEqual(out["counts"]["library_primary_records"], str(PRIMARY))
                self.assertEqual(out["aligner"]["low_as_mode"], "off")
                self.assertIsNone(out["ids"])

    def test_off_is_byte_identical_to_plain_flag4_extraction(self):
        """Off, the FASTQ equals `samtools view -F 2304 -f 4 | samtools fastq -N` -- the pre-rescue path."""
        out = self.run_module("auto", (PG_STAR,))
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        (tmp / "in.sam").write_text(sam_text((PG_STAR,)))
        tool = str(SAMTOOLS_BIN / "samtools")
        view = subprocess.run([tool, "view", "-F", "2304", "-f", "4", "-u", str(tmp / "in.sam")],
                              capture_output=True, check=True).stdout
        fq = subprocess.run([tool, "fastq", "-N", "-"], input=view, capture_output=True, check=True).stdout
        lines = fq.decode().splitlines()
        plain = sorted((lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4))
        self.assertEqual(out["reads"], plain)

    def test_on_overrides_a_star_header_and_off_a_bwa_T0_header(self):
        on = self.run_module("on", (PG_STAR,))
        self.assertEqual(on["reads"], expected_reads(("flag4", "low")))
        self.assertEqual(on["aligner"]["low_as_mode_reason"], "requested_on")
        off = self.run_module("off", (PG_GDC,))
        self.assertEqual(off["reads"], expected_reads(("flag4",)))
        self.assertEqual(off["aligner"]["low_as_mode_reason"], "requested_off")

    def test_cram_with_reference(self):
        out = self.run_module("auto", (PG_GDC,), kind="cram")
        self.assertEqual(out["reads"], expected_reads(("flag4", "low")))
        self.assertEqual(out["counts"]["low_as_reads_extracted"], "3")
        self.assertEqual(out["counts"]["library_primary_records"], str(PRIMARY))
        off = self.run_module("auto", (PG_STAR,), kind="cram")
        self.assertEqual(off["reads"], expected_reads(("flag4",)))

    def test_provenance_never_records_the_command_line(self):
        out = self.run_module("auto", (PG_GDC, PG_SAMTOOLS))
        self.assertEqual((out["aligner"]["upstream_aligner"], out["aligner"]["upstream_aligner_version"],
                          out["aligner"]["upstream_bwa_min_score"]), ("bwa", "0.7.15-r1140", "0"))
        self.assertNotIn("SECRET", "".join(out["aligner"].values()) + out["stderr"])

    def test_no_scratch_files_left_behind(self):
        for mode in ("on", "off"):
            with self.subTest(mode=mode):
                left = self.run_module(mode, (PG_GDC,))["leftovers"]
                allowed = {".command.sh", "S.UNMAPPED.fastq.gz", "S.extract_counts.tsv",
                           "S.upstream_aligner.tsv", "S.low_as_read_ids.txt.gz"}
                self.assertEqual(set(left) - allowed, set())


# ------------------------------------------------------------------ the auto decision


class AutoDecision(unittest.TestCase):
    def decide(self, *pg, mode="auto"):
        return LSR.decide(header(*pg), mode, "S")

    def test_bwa_mem_below_30_turns_auto_on(self):
        for pg, score in ((PG_GDC, "0"), (PG_BWA_T19, "19"), (PG_BWA_MEM2, "0"), (PG_SENTIEON, "0")):
            with self.subTest(pg=pg):
                row = self.decide(pg)
                self.assertEqual((row["low_as_mode"], row["low_as_mode_reason"], row["upstream_bwa_min_score"]),
                                 ("on", "auto_bwa_mem_T_below_30", score))

    def test_bwa_mem_at_default_or_30_leaves_auto_off(self):
        self.assertEqual(self.decide(PG_BWA_DEFAULT)["low_as_mode_reason"], "auto_bwa_mem_T_default")
        self.assertEqual(self.decide(PG_BWA_DEFAULT)["upstream_bwa_min_score"], "default")
        self.assertEqual(self.decide(PG_BWA_T30)["low_as_mode_reason"], "auto_bwa_mem_T_30_or_more")
        for pg in (PG_BWA_DEFAULT, PG_BWA_T30):
            self.assertEqual(self.decide(pg)["low_as_mode"], "off")

    def test_star_minimap2_bwa_aln_and_no_pg_leave_auto_off(self):
        cases = {(PG_STAR,): ("auto_other_aligner", "STAR", "2.7.5c"),
                 (PG_MM2,): ("auto_other_aligner", "minimap2", "2.28-r1209"),
                 (PG_ALN,): ("auto_bwa_not_mem", "bwa", "0.7.12"),
                 (PG_NO_CL,): ("auto_bwa_not_mem", "bwa", "0.7.17"),
                 (): ("auto_no_pg", "NA", "NA"),
                 (PG_UNKNOWN,): ("auto_unknown_aligner", "mytool", "1")}
        for pg, (reason, aligner, version) in cases.items():
            with self.subTest(pg=pg):
                row = self.decide(*pg)
                self.assertEqual(row["low_as_mode"], "off")
                self.assertEqual((row["low_as_mode_reason"], row["upstream_aligner"],
                                  row["upstream_aligner_version"], row["upstream_bwa_min_score"]),
                                 (reason, aligner, version, "NA"))

    def test_on_and_off_override_detection(self):
        self.assertEqual(self.decide(PG_STAR, mode="on")["low_as_mode"], "on")
        self.assertEqual(self.decide(PG_STAR, mode="on")["upstream_aligner"], "STAR")   # still recorded
        self.assertEqual(self.decide(PG_GDC, mode="off")["low_as_mode"], "off")
        self.assertEqual(self.decide(PG_GDC, mode="off")["upstream_bwa_min_score"], "0")
        with self.assertRaises(ValueError):
            self.decide(PG_GDC, mode="maybe")

    def test_minus_t_forms(self):
        self.assertEqual(LSR.bwa_mem_min_score("bwa mem -T0 r.fa"), "0")
        self.assertEqual(LSR.bwa_mem_min_score("bwa mem -t 4 -T 25 r.fa"), "25")
        self.assertEqual(LSR.bwa_mem_min_score("bwa mem -T r.fa"), "unparsed")
        self.assertEqual(LSR.bwa_mem_min_score("bwa mem r.fa"), "default")
        self.assertIsNone(LSR.bwa_mem_min_score("bwa aln r.fa a.fq"))
        self.assertEqual(self.decide("@PG\tID:bwa\tPN:bwa\tCL:bwa mem -T r.fa")["low_as_mode_reason"],
                         "auto_bwa_mem_T_unparsed")

    def test_command_line_is_never_recorded(self):
        for pg in (PG_GDC, PG_BWA_DEFAULT, PG_STAR):
            row = self.decide(pg, PG_SAMTOOLS)
            self.assertNotIn("SECRET", "".join(map(str, row.values())))
            self.assertNotIn("SECRET", LSR.explain(row))

    def test_every_reason_maps_to_on_or_off(self):
        self.assertEqual(set(LSR.REASONS.values()), {"on", "off"})
        self.assertEqual(QC.LOW_AS_MODES, {"on", "off"})


# ------------------------------------------------------------------ survival counting


class Survival(unittest.TestCase):
    def test_counts_rescued_names_after_depletion_and_classification(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with gzip.open(d / "ids.txt.gz", "wt") as fh:
                fh.write("r1/1\nr1/2\nr2\nr3/1\n")
            with gzip.open(d / "dep.fastq.gz", "wt") as fh:
                for name in ("r1/1", "r1/2", "r2", "flag4_a/1", "flag4_b"):
                    fh.write(f"@{name} comment\nACGT\n+\nIIII\n")
            (d / "cls.fasta").write_text(">r1/2 kraken:taxid|3001\nACGT\n>flag4_a/1\nACGT\n")
            row = LSR.survival(d / "ids.txt.gz", d / "dep.fastq.gz", d / "cls.fasta", "S")
            self.assertEqual(row, {"sample": "S", "low_as_reads_after_depletion": 3, "low_as_reads_classified": 1})
            self.assertEqual(LSR.survival(d / "ids.txt.gz", d / "dep.fastq.gz", None, "S")["low_as_reads_classified"], "NA")

    def test_empty_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with gzip.open(d / "ids.txt.gz", "wt"):
                pass
            with gzip.open(d / "dep.fastq.gz", "wt"):
                pass
            (d / "cls.fasta").write_text("")
            row = LSR.survival(d / "ids.txt.gz", d / "dep.fastq.gz", d / "cls.fasta", "S")
            self.assertEqual((row["low_as_reads_after_depletion"], row["low_as_reads_classified"]), (0, 0))


# ------------------------------------------------------------------ QC summary and run record


def write(path, text):
    path.write_text(text)


class QcSummaryColumns(unittest.TestCase):
    def build(self, tmp, extract, aligner=None, survival=None, kraken_root=None):
        d = Path(tmp)
        write(d / "B.extract_counts.tsv",
              "sample\tlibrary_primary_records\textracted_unmapped_records\tlow_as_reads_extracted\tlow_as_no_as_records\n"
              f"B\t{extract}\n")
        if aligner:
            write(d / "B.upstream_aligner.tsv", "\t".join(LSR.DECISION_FIELDS) + "\n" + "B\t" + aligner + "\n")
        if survival:
            write(d / "B.low_as_survival.tsv", "\t".join(LSR.SURVIVAL_FIELDS) + "\n" + "B\t" + survival + "\n")
        if kraken_root is not None:
            write(d / "B.krakenuniq.report.txt",
                  f"50\t10\t10\t0\t0\tNA\t0\tno rank\tunclassified\n50\t{kraken_root}\t0\t0\t0\t0\t1\tno rank\troot\n")
        return {r["sample"]: r for r in QC.build(d, 100)}["B"]

    def test_columns_present(self):
        for column in ("upstream_aligner", "upstream_aligner_version", "upstream_bwa_min_score", "low_as_mode",
                       "low_as_mode_reason", "low_as_reads_extracted", "low_as_no_as_records",
                       "low_as_reads_after_depletion", "low_as_reads_classified", "low_as_share_of_bracken"):
            self.assertIn(column, QC.COLUMNS)

    def test_on_sample_reports_each_stage_and_share(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = self.build(tmp, "100\t40\t8\t1", "bwa\t0.7.15-r1140\t0\ton\tauto_bwa_mem_T_below_30",
                             "6\t4", kraken_root=40)
            self.assertEqual((row["extracted_unmapped_records"], row["low_as_reads_extracted"],
                              row["low_as_reads_after_depletion"], row["low_as_reads_classified"]),
                             ("40", "8", "6", "4"))
            self.assertEqual(row["low_as_share_of_bracken"], "10.00")   # 4 of 40 classified
            self.assertEqual(row["pct_library_unmapped"], "40.00")       # FLAG 4 only, unchanged meaning
            self.assertEqual(row["upstream_aligner"], "bwa")

    def test_off_sample_reads_zero_not_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            row = self.build(tmp, "100\t40\t0\tNA", "STAR\t2.7.5c\tNA\toff\tauto_other_aligner", kraken_root=20)
            self.assertEqual((row["low_as_reads_extracted"], row["low_as_reads_after_depletion"],
                              row["low_as_reads_classified"], row["low_as_share_of_bracken"]),
                             ("0", 0, 0, "0.00"))
            self.assertEqual(row["low_as_no_as_records"], "NA")

    def test_pre_rescue_fragment_still_parses(self):
        with tempfile.TemporaryDirectory() as tmp:
            write(Path(tmp) / "B.extract_counts.tsv",
                  "sample\tlibrary_primary_records\textracted_unmapped_records\nB\t31\t28\n")
            row = {r["sample"]: r for r in QC.build(Path(tmp), 100)}["B"]
            self.assertEqual((row["extracted_unmapped_records"], row["low_as_reads_extracted"], row["low_as_mode"]),
                             ("28", "NA", "NA"))

    def test_inconsistent_counts_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "exceeds"):
                self.build(tmp, "100\t40\t8\t1", "bwa\t0.7.15\t0\ton\tauto_bwa_mem_T_below_30", "9\t4")
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "low_as_mode off"):
                self.build(tmp, "100\t40\t8\tNA", "STAR\t2.7.5c\tNA\toff\tauto_other_aligner")
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "on or off"):
                self.build(tmp, "100\t40\t8\tNA", "STAR\t2.7.5c\tNA\tyes\tauto_other_aligner")


class RunRecord(unittest.TestCase):
    def test_cohort_aggregate_is_counts_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            qc = Path(tmp)
            cols = ["sample", "upstream_aligner", "upstream_aligner_version", "upstream_bwa_min_score",
                    "low_as_mode", "low_as_mode_reason", "low_as_reads_extracted",
                    "low_as_reads_after_depletion", "low_as_reads_classified"]
            rows = [["TCGA-SECRET-1", "bwa", "0.7.15", "0", "on", "auto_bwa_mem_T_below_30", "8", "6", "4"],
                    ["TCGA-SECRET-2", "bwa", "0.7.15", "0", "on", "auto_bwa_mem_T_below_30", "2", "1", "1"],
                    ["TCGA-SECRET-3", "STAR", "2.7.5c", "NA", "off", "auto_other_aligner", "0", "0", "0"],
                    ["FASTQ-SECRET", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]]
            (qc / "cmp.qc.summary.tsv").write_text(
                "\n".join("\t".join(r) for r in [cols] + rows) + "\n")
            summary = FP.low_score_rescue_summary({"qc_summary_dir": str(qc)})
            self.assertEqual(summary["inputs"], 3)
            self.assertEqual(summary["by_mode"], {"off": 1, "on": 2})
            self.assertEqual(summary["by_reason"], {"auto_bwa_mem_T_below_30": 2, "auto_other_aligner": 1})
            self.assertEqual(summary["rescued_reads_total"]["low_as_reads_extracted"], 10)
            self.assertEqual(summary["rescued_reads_total"]["low_as_reads_classified"], 5)
            self.assertNotIn("SECRET", repr(summary))
            report = FP.render({"status": "completed", "outputs": {"low_score_rescue": summary}},
                               [], {}, Path("/t"), "s")
            self.assertIn("LOW-SCORE RESCUE", report)
            self.assertNotIn("SECRET", report)
            self.assertIsNone(FP.low_score_rescue_summary({"qc_summary_dir": str(qc / "missing")}))


# ------------------------------------------------------------------ wiring and docs


class Wiring(unittest.TestCase):
    def test_param_default_and_strict_check(self):
        self.assertIn('params.extract_low_score_mapped = "auto"', MAIN)
        self.assertIn("--extract_low_score_mapped must be auto, on or off", MAIN)
        self.assertIn("instanceof CharSequence", MAIN)   # a bare flag (Boolean true) is refused

    def test_threshold_is_one_number(self):
        """AS < 30 in the module's samtools expression and awk equals the script's constant."""
        self.assertEqual(LSR.LOW_AS_THRESHOLD, 30)
        self.assertIn(f"[AS] < {LSR.LOW_AS_THRESHOLD}'", MODULE)
        self.assertIn(f"-v min_as={LSR.LOW_AS_THRESHOLD}", MODULE)

    def test_decided_from_the_header_before_the_decode(self):
        self.assertIn('--mode "${params.extract_low_score_mapped}"', MODULE)
        self.assertLess(MODULE.index("low_score_rescue.py\" decide"), MODULE.index("mkfifo primary.fifo"))

    def test_one_decode_per_branch_with_cram_reference(self):
        on = MODULE[MODULE.index('if [[ "\\$LOW_AS" == "on" ]]'):MODULE.index("\n    else\n")]
        off = MODULE[MODULE.index("\n    else\n"):MODULE.index("\n    fi\n    printf")]
        for branch in (on, off):
            self.assertEqual(branch.count('"${alignment}"'), 1)
            self.assertIn("-F 2304", branch)
            self.assertIn('"\\${REFERENCE_ARGS[@]}"', branch)

    def test_off_branch_is_the_pre_rescue_pipeline(self):
        """Off, the extraction pipeline is the one CMPipeline ran before the rescue, verbatim."""
        before = [
            'samtools view -c unmapped.fifo > unmapped.count & P_UNMAPPED=\\$!',
            'if ! samtools view -@ "${task.cpus}" -F 2304 -u \\',
            '"\\${REFERENCE_ARGS[@]}" "${alignment}" |',
            'tee primary.fifo |',
            'samtools view -f 4 -u - |',
            'tee unmapped.fifo |',
            'samtools fastq -@ "${task.cpus}" -N - |',
            'bgzip -@ "${task.cpus}" -c > "\\$READS"; then',
        ]
        off = MODULE[MODULE.index("\n    else\n"):MODULE.index("\n    fi\n    printf")]
        lines = [line.strip() for line in off.splitlines() if line.strip()]
        start = lines.index(before[0])
        self.assertEqual(lines[start:start + len(before)], before)

    def test_survival_is_wired(self):
        self.assertIn("include { lowAsSurvival } from './Modules/low_as_survival.nf'", MAIN)
        self.assertIn("extractReads.out.low_as_ids.join(MAPPED_READS_MULTI.PAN)", MAIN)
        self.assertIn("qcParts << lowAsSurvival.out.counts", MAIN)
        self.assertIn("qcParts << extractReads.out.aligner", MAIN)
        self.assertIn('"${params.scripts}/low_score_rescue.py"', MAIN)   # in depsOf: -resume invalidates
        self.assertIn("low_score_rescue.py\" survival", SURVIVAL_MODULE)

    def test_id_list_is_not_published_by_default(self):
        self.assertIn("fn.endsWith('.low_as_read_ids.txt.gz')", MODULE)

    def test_documented(self):
        for path in ("docs/parameters.md", "docs/running/classification.md", "docs/CHANGELOG.md", "docs/outputs.md"):
            with self.subTest(path=path):
                text = (ROOT / path).read_text()
                self.assertIn("extract_low_score_mapped", text)
        self.assertIn("-T 0", (ROOT / "docs/CHANGELOG.md").read_text())


if __name__ == "__main__":
    unittest.main()
