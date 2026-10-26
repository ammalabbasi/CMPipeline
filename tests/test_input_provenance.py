#!/usr/bin/env python3
"""Input provenance from BAM/CRAM headers (audit X13, Ludmil E3; pksProfiler 70cb629).

Synthetic headers only: invented contig sets, invented read groups and sample names (SECRET_*),
which must never reach any output. The real-run path (extractReads writes the provenance row from
the header it already reads) is exercised in tests/test_low_score_rescue.py's extraction tests.
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HP = load("header_provenance", ROOT / "scripts" / "header_provenance.py")
QC = load("build_qc_summary_hp", ROOT / "scripts" / "build_qc_summary.py")
FP = load("finalize_provenance_hp", ROOT / "scripts" / "finalize_provenance.py")
LSR = load("low_score_rescue_hp", ROOT / "scripts" / "low_score_rescue.py")

GRCH38_CHR1 = 248956422


def sq(name, length, m5=None):
    return f"@SQ\tSN:{name}\tLN:{length}" + (f"\tM5:{m5}" if m5 else "")


# A GDC GRCh38.d1.vd1-like contig set, in miniature: primary, ALT, random, decoy, EBV, HPV, HLA.
GDC_SQ = [sq("chr1", GRCH38_CHR1), sq("chr2", 242193529), sq("chrM", 16569),
          sq("chr1_KI270706v1_random", 175055), sq("chr1_KI270762v1_alt", 354444),
          sq("chrUn_JTFH01000001v1_decoy", 25139), sq("chrUn_JTFH01000002v1_decoy", 18532),
          sq("chrEBV", 171823), sq("HPV16", 7906), sq("HPV-mSD2", 7300), sq("HCV-1", 9646),
          sq("HLA-A*01:01:01:01", 3503)]
PLAIN_SQ = [sq("chr1", GRCH38_CHR1), sq("chr2", 242193529), sq("chrM", 16569)]
GRCH37_SQ = [sq("1", 249250621), sq("2", 243199373), sq("MT", 16569)]

PG_GDC = ("@PG\tID:bwa\tPN:bwa\tVN:0.7.15-r1140\tCL:bwa mem -t 36 -T 0 -R "
          "@RG\\tID:SECRET_RG\\tSM:SECRET_SAMPLE\\tPU:SECRET_FLOWCELL.1 /SECRET/GRCh38.d1.vd1.fa "
          "/SECRET/SECRET_1.fq /SECRET/SECRET_2.fq")
PG_SORT = "@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.3\tCL:samtools sort -o /SECRET/SECRET.bam"
PG_MARKDUP = ("@PG\tID:MarkDuplicates\tPN:MarkDuplicates\tPP:samtools\tVN:2.4.1\t"
              "CL:picard MarkDuplicates INPUT=/SECRET/x.bam OUTPUT=/SECRET/y.bam")
RG = "@RG\tID:SECRET_RG\tSM:SECRET_SAMPLE\tLB:SECRET_LIB\tPU:SECRET_FLOWCELL.1"


def header(sqs, *pg, rg=True):
    return "\n".join(["@HD\tVN:1.6\tSO:coordinate"] + sqs + ([RG] if rg else []) + list(pg)) + "\n"


class Fingerprint(unittest.TestCase):
    def test_same_reference_same_fingerprint_whatever_the_pg_and_rg(self):
        a = HP.parse(header(GDC_SQ, PG_GDC, PG_SORT), "A")
        b = HP.parse(header(GDC_SQ, PG_GDC, rg=False), "B")
        self.assertEqual(a["input_ref_fingerprint"], b["input_ref_fingerprint"])

    def test_different_builds_and_contig_sets_differ(self):
        fps = {HP.parse(header(s), "x")["input_ref_fingerprint"] for s in (GDC_SQ, PLAIN_SQ, GRCH37_SQ)}
        self.assertEqual(len(fps), 3)
        # Order is part of the reference identity (it is how the BAM's tid map works).
        self.assertNotEqual(HP.parse(header(PLAIN_SQ), "x")["input_ref_fingerprint"],
                            HP.parse(header(PLAIN_SQ[::-1]), "x")["input_ref_fingerprint"])

    def test_m5_fingerprint_needs_every_m5_and_ignores_names(self):
        self.assertEqual(HP.parse(header(PLAIN_SQ), "x")["input_ref_m5_fingerprint"], "NA")
        ucsc = [sq("chr1", 10, "a" * 32), sq("chr2", 20, "b" * 32)]
        ensembl = [sq("1", 10, "a" * 32), sq("2", 20, "b" * 32)]
        u, e = HP.parse(header(ucsc), "x"), HP.parse(header(ensembl), "x")
        self.assertEqual(u["input_ref_m5_fingerprint"], e["input_ref_m5_fingerprint"])
        self.assertNotEqual(u["input_ref_fingerprint"], e["input_ref_fingerprint"])

    def test_build_and_contig_classes(self):
        row = HP.parse(header(GDC_SQ, PG_GDC), "S")
        self.assertEqual(row["input_ref_build"], "GRCh38")
        self.assertEqual(row["input_ref_contigs"], len(GDC_SQ))
        self.assertEqual((row["input_ref_alt_contigs"], row["input_ref_decoy_contigs"],
                          row["input_ref_hla_contigs"], row["input_ref_viral_contigs"]), (1, 2, 1, 4))
        self.assertEqual(row["input_ref_has_ebv"], "true")
        plain = HP.parse(header(PLAIN_SQ), "S")
        self.assertEqual((plain["input_ref_decoy_contigs"], plain["input_ref_viral_contigs"],
                          plain["input_ref_has_ebv"]), (0, 0, "false"))
        self.assertEqual(HP.parse(header(GRCH37_SQ), "S")["input_ref_build"], "GRCh37")
        self.assertEqual(HP.parse(header([sq("contig_a", 5)]), "S")["input_ref_build"], "no_chr1")

    def test_no_sq_is_a_named_error(self):
        with self.assertRaisesRegex(ValueError, "no @SQ"):
            HP.parse(header([], PG_GDC), "S")


class AlignerChain(unittest.TestCase):
    def test_chain_in_pp_order_program_version_only(self):
        row = HP.parse(header(GDC_SQ, PG_MARKDUP, PG_SORT, PG_GDC), "S")   # listed out of order
        self.assertEqual(row["input_pg_chain"], "bwa:0.7.15-r1140 > samtools:1.3 > MarkDuplicates:2.4.1")
        self.assertEqual(HP.parse(header(GDC_SQ), "S")["input_pg_chain"], "NA")

    def test_bwa_options_reduced_to_letters_and_numbers(self):
        self.assertEqual(HP.bwa_mem_options(PG_GDC.split("CL:", 1)[1]), "-T 0")
        self.assertEqual(HP.bwa_mem_options(
            "bwa mem -K 100000000 -Y -M -T30 -O 6,6 -x ont2d -o /SECRET/out.sam ref a.fq"),
            "-K 100000000 -M -O 6,6 -T 30 -Y")
        self.assertEqual(HP.bwa_mem_options("bwa mem ref.fa a.fq"), "default")
        self.assertEqual(HP.bwa_mem_options("bwa-mem2 mem -t 8 -T 0 ref a b"), "-T 0")
        self.assertEqual(HP.bwa_mem_options("bwa mem -k SECRETVALUE ref a"), "-k ?")
        self.assertIsNone(HP.bwa_mem_options("STAR --runThreadN 8"))

    def test_nothing_identifying_reaches_the_row(self):
        row = HP.parse(header(GDC_SQ, PG_GDC, PG_SORT, PG_MARKDUP), "S")
        text = "\t".join(str(v) for v in row.values())
        for secret in ("SECRET", "/", "@RG", "ID:", "SM:", "PU:", "INPUT="):
            self.assertNotIn(secret, text)


class UnmappedFilter(unittest.TestCase):
    CASES = [
        ("samtools view -b -F 4 -o out.bam in.bam", {"samtools_view_exclude_unmapped"}),
        ("samtools view -b -F0x904 in.bam", {"samtools_view_exclude_unmapped"}),
        ("samtools view -b --exclude-flags UNMAP,SECONDARY in.bam", {"samtools_view_exclude_unmapped"}),
        ("samtools view -b -F 0x900 in.bam", set()),                 # secondary/supplementary only
        ("samtools view -b -F 1024 in.bam", set()),                  # duplicates only
        ("samtools view -b -f 2 in.bam", {"samtools_view_require_proper_pair"}),
        ("samtools view -b -q 20 in.bam", {"samtools_view_min_mapq"}),
        ("samtools view -b -q 0 in.bam", set()),
        ("samtools view -b -L targets.bed in.bam", {"samtools_view_region"}),
        ("samtools view -b in.bam chr1:1-100", {"samtools_view_region"}),
        ("samtools view -b --input-fmt-option decode_md=0 -o out.bam in.bam", set()),
        ("samtools view -e '!flag.unmap' in.bam", {"samtools_view_exclude_unmapped"}),
        ("sambamba view -F 'not unmapped' in.bam", {"sambamba_view_filter"}),
        ("samtools sort -o x.bam in.bam", set()),
        ("bwa mem -T 0 ref a b", set()),
    ]

    def test_cases(self):
        for cl, expected in self.CASES:
            with self.subTest(cl=cl):
                self.assertEqual(HP.unmapped_filters(cl), expected)
                self.assertTrue(expected <= set(HP.FILTER_REASONS))

    def test_row_value(self):
        pg = "@PG\tID:samtools.1\tPN:samtools\tPP:bwa\tVN:1.9\tCL:samtools view -b -F 4 -o /SECRET/x.bam"
        self.assertEqual(HP.parse(header(GDC_SQ, PG_GDC, pg), "S")["input_unmapped_filter"],
                         "samtools_view_exclude_unmapped")
        self.assertEqual(HP.parse(header(GDC_SQ, PG_GDC, PG_SORT), "S")["input_unmapped_filter"], "none")


class DecideWritesProvenanceFromTheSameHeader(unittest.TestCase):
    def test_decide_cli_writes_both_rows(self):
        import subprocess
        import sys
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "h.sam").write_text(header(GDC_SQ, PG_GDC, PG_SORT))
            r = subprocess.run([sys.executable, str(ROOT / "scripts" / "low_score_rescue.py"), "decide",
                                "--header", str(d / "h.sam"), "--sample", "S", "--mode", "auto",
                                "--output", str(d / "a.tsv"), "--provenance", str(d / "p.tsv")],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(r.stdout.strip(), "on")
            with open(d / "p.tsv") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(list(rows[0]), HP.PROVENANCE_FIELDS)
            self.assertEqual(rows[0]["input_bwa_options"], "-T 0")
            self.assertNotIn("SECRET", (d / "p.tsv").read_text() + r.stderr)


class QcAndRunReport(unittest.TestCase):
    def build_cohort(self, d: Path):
        inputs = d / "in"
        inputs.mkdir()
        for i, (sqs, pgs) in enumerate([(GDC_SQ, (PG_GDC, PG_SORT))] * 3 + [(PLAIN_SQ, (PG_SORT,))]):
            name = f"SECRET_S{i}"
            HP.write_row(inputs / f"{name}.input_provenance.tsv", HP.parse(header(sqs, *pgs), name))
            LSR.write_row(inputs / f"{name}.upstream_aligner.tsv", LSR.DECISION_FIELDS,
                          LSR.decide(header(sqs, *pgs), "auto", name))
        rows = QC.build(inputs, 100000)
        out = d / "QC"
        out.mkdir()
        with open(out / "cmp.qc.summary.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=QC.COLUMNS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        return rows, out

    def test_qc_columns_and_cohort_summary_counts_only(self):
        with tempfile.TemporaryDirectory() as d:
            rows, out = self.build_cohort(Path(d))
            for column in HP.PROVENANCE_FIELDS[1:]:
                self.assertIn(column, QC.COLUMNS)
            self.assertEqual({r["input_ref_build"] for r in rows}, {"GRCh38"})
            summary = FP.input_alignments_summary({"qc_summary_dir": str(out)})
            self.assertEqual(summary["inputs"], 4)
            self.assertEqual(len(summary["references"]), 2)
            self.assertEqual(summary["aligners"], {"bwa 0.7.15-r1140 -T 0": 3, "NA NA": 1})
            lines = "\n".join(FP.render_input_alignments(summary))
            self.assertIn("4 input(s): 2 reference fingerprint(s); 2 aligner configuration(s)", lines)
            self.assertIn("more than one input reference", lines)
            self.assertNotIn("SECRET", json.dumps(summary) + lines)

    def test_fastq_only_run_has_no_section(self):
        self.assertIsNone(FP.input_alignments_summary({"qc_summary_dir": "/nonexistent"}))
        self.assertEqual(FP.render_input_alignments(None), [])


class Wiring(unittest.TestCase):
    def test_module_and_main(self):
        module = (ROOT / "Modules" / "extract_reads.nf").read_text()
        main = (ROOT / "main.nf").read_text()
        self.assertIn('--provenance "${sampleID}.input_provenance.tsv"', module)
        self.assertIn('path("${sampleID}.input_provenance.tsv"), emit: provenance', module)
        # One header read from the alignment per sample.
        self.assertEqual(module.count(" -H "), 1)
        self.assertIn("qcParts << extractReads.out.provenance", main)
        self.assertIn('"${params.scripts}/header_provenance.py"', main)

    def test_column_count_documented(self):
        n = len(QC.COLUMNS)
        for path in ("docs/outputs.md", "README.md"):
            with self.subTest(path=path):
                self.assertIn(f"{n} columns", (ROOT / path).read_text())

    def test_documented(self):
        for path in ("docs/outputs.md", "docs/CHANGELOG.md"):
            with self.subTest(path=path):
                text = (ROOT / path).read_text()
                self.assertIn("input_ref_fingerprint", text)
                self.assertIn("input_unmapped_filter", text)


if __name__ == "__main__":
    unittest.main()
