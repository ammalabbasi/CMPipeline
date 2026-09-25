#!/usr/bin/env python3
"""main.nf metadataColumns(): the launch-time metadata header parser.

Expected values come from the original quote-toggling scan (the pre-rewrite for-loop version).
The rewrite must keep those results and must not use loops, which the Nextflow 26 strict parser
rejects. The executed check needs nextflow on PATH (it runs the function body taken from main.nf).
"""
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CASES = [
    ["sampleid\tsample_type\tshipment_batch\tage_diag\tsex", ["sampleid", "sample_type", "shipment_batch", "age_diag", "sex"]],
    ["patient,\"Sample Type_x\",batch", ["patient", "Sample Type_x", "batch"]],
    ["\"a,b\",c,\"d\"\"e\"", ["a,b", "c", "d\"e"]],
    ["x\t\"y z\"\tw", ["x", "y z", "w"]],
    ["a\"b,c\"d,e", ["ab,cd", "e"]],
    ["a\"b\"\"c\"d,e", ["ab\"cd", "e"]],
    ["\"\t\"aa,,\"\"", ["\taa,,"]],
    [",\"\"\"\t\"aab", [",\"\taab"]],
    ["a\",a\"\"\"\ta", ["a,a\"", "a"]],
    ["\"a\ta\",\"a\"", ["a\ta,a"]],
    ["a,\"b\"c\"d\",e", ["a", "bcd", "e"]],
    ["a,,b", "ERR_BLANKDUP"],
    ["a,a", "ERR_BLANKDUP"],
    ["\"\"", "ERR_BLANKDUP"],
    ["\"a,b", "ERR_UNCLOSED"],
]


def function_source():
    text = (ROOT / "main.nf").read_text()
    start = text.index("def metadataColumns(String header) {")
    return text[start:text.index("\n}\n", start) + 3]


class MetadataHeaderParser(unittest.TestCase):
    def test_no_loops_for_the_strict_parser(self):
        body = function_source()
        self.assertNotRegex(body, r"\bfor\s*\(|\bwhile\s*\(")

    @unittest.skipUnless(shutil.which("nextflow"), "nextflow not on PATH")
    def test_matches_the_original_scan(self):
        body = (function_source()
                .replace('error "Unclosed quote in metadata header"', 'return "ERR_UNCLOSED"')
                .replace('error "Blank or duplicate metadata header"', 'return "ERR_BLANKDUP"'))
        self.assertNotIn("error ", body)
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "cases.json").write_text(json.dumps(CASES))
            (d / "t.nf").write_text(body + '''
workflow {
    def cases = new groovy.json.JsonSlurper().parseText(file("cases.json").text)
    def bad = cases.findAll { c -> metadataColumns(c[0]) != c[1] }
    println "MISMATCHES=${bad.size()} ${groovy.json.JsonOutput.toJson(bad)}"
}
''')
            r = subprocess.run(["nextflow", "-q", "run", "t.nf"], cwd=d, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("MISMATCHES=0 []", r.stdout, r.stdout)


if __name__ == "__main__":
    unittest.main()
