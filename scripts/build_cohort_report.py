#!/usr/bin/env python3
"""Cohort report: one self-contained HTML page summarising a CMPipeline run.

Reads the QC summary table (scripts/build_qc_summary.py) and, when present, the run
provenance record. It deliberately does not re-walk the results tree: the summary table
is the single source, so the report and the table can never disagree.

The page embeds everything -- no network, no CDN, no JavaScript libraries -- so it opens
from a cluster filesystem or an email attachment years from now.

Counts and rates only. No sequence and no clinical value is read or written. The page
DOES carry per-sample identifiers, so treat it as you would any sample-level table.

Colour: the depletion stages are an ordinal ramp (ordered funnel stages), rendered as a
single blue hue light->dark, with the retained segment in a distinct hue because it is
the outcome rather than another stage. Light and dark steps were each chosen against
their own surface and validated for adjacent-pair separation, in normal vision and under
protanopia and deuteranopia simulation; worst adjacent pair light dE 29.0 normal /
28.7 CVD, dark 19.2 / 19.2, against floors of 15 and 8. Two light swatches sit under
3:1 contrast, so every segment is directly labelled and the same numbers appear in a
table below the chart -- identity is never carried by colour alone.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
from datetime import datetime, timezone
from pathlib import Path

NA = {"", "NA", None}

STAGES = [
    ("pct_removed_hg38", "Removed: GRCh38", "--seg-1"),
    ("pct_removed_t2t_phix", "Removed: T2T + PhiX", "--seg-2"),
    ("pct_removed_pangenome", "Removed: pangenome", "--seg-3"),
    ("__retained", "Retained (microbial candidate)", "--seg-4"),
]


def num(value, default=None):
    if value in NA:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def fmt_int(value):
    n = num(value)
    return "—" if n is None else f"{int(n):,}"


def fmt_pct(value, digits=1):
    n = num(value)
    return "—" if n is None else f"{n:.{digits}f}%"


def median(values):
    values = sorted(v for v in values if v is not None)
    if not values:
        return None
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2


# ------------------------------------------------------------------------ sections


def stat_tiles(rows):
    total = len(rows)
    raw = sum(num(r.get("reads_raw"), 0) for r in rows)
    retained = median([num(r.get("pct_retained_overall")) for r in rows])
    gate = [r.get("cleared_consensus_gate") for r in rows]
    cleared = sum(1 for g in gate if g == "true")
    known = sum(1 for g in gate if g in ("true", "false"))
    # BAM/CRAM samples enter QC with only their unmapped reads (audit C07), so "reads entering
    # QC" is not the library size. When extraction counts exist, show the library too.
    library = sum(num(r.get("library_primary_records"), 0) for r in rows)
    lib_retained = median([num(r.get("pct_retained_of_library")) for r in rows])

    tiles = [
        ("Samples", f"{total:,}", "in this cohort"),
        ("Reads entering QC", f"{raw:,.0f}" if raw else "—",
         "all samples; for BAM/CRAM, only the unmapped reads extracted"),
        ("Median retention",
         "—" if retained is None else f"{retained:.1f}%",
         "of reads entering QC that survive host depletion"),
    ]
    if library:
        tiles.append(("Library records (BAM/CRAM)", f"{library:,.0f}",
                      "primary records in the input alignments"))
        tiles.append(("Median retention of library",
                      "—" if lib_retained is None else f"{lib_retained:.3f}%",
                      "BAM/CRAM samples: host-depleted reads / primary records"))
    tiles += [
        ("Cleared the consensus gate",
         "—" if not known else f"{cleared} of {known}",
         "eligible for MetaPhlAn"),
    ]
    cells = "".join(
        f'<div class="tile"><div class="tile-label">{html.escape(label)}</div>'
        f'<div class="tile-value">{html.escape(value)}</div>'
        f'<div class="tile-note">{html.escape(note)}</div></div>'
        for label, value, note in tiles
    )
    return f'<section class="tiles">{cells}</section>'


def retention_chart(rows):
    """Stacked horizontal bars: where each sample's reads went.

    Every segment is labelled where it is wide enough, and a table below carries the
    same numbers, so the chart is never the only way to read the data.
    """
    usable = [r for r in rows if num(r.get("pct_retained_overall")) is not None]
    if not usable:
        return ('<section><h2>Read retention</h2><p class="empty">No depletion counts in '
                'the summary table — host depletion did not run, or predates the '
                '<code>depletion_counts.tsv</code> output.</p></section>')

    row_h, gap, label_w, bar_w = 26, 6, 190, 560
    height = len(usable) * (row_h + gap) + 34
    parts = []

    for i, row in enumerate(usable):
        y = i * (row_h + gap) + 28
        sample = row.get("sample", "?")
        retained = num(row.get("pct_retained_overall"), 0.0)
        segs = []
        for key, label, token in STAGES:
            value = retained if key == "__retained" else num(row.get(key), 0.0)
            segs.append((label, max(0.0, value or 0.0), token))
        total = sum(v for _, v, _ in segs) or 100.0

        parts.append(
            f'<text class="cat" x="{label_w - 10}" y="{y + row_h / 2 + 4}" '
            f'text-anchor="end">{html.escape(sample)}</text>'
        )
        x = label_w
        for label, value, token in segs:
            w = (value / total) * bar_w
            if w <= 0:
                continue
            # 2px surface gap between adjacent fills
            draw_w = max(0.0, w - 2)
            parts.append(
                f'<rect class="seg" x="{x:.1f}" y="{y}" width="{draw_w:.1f}" height="{row_h}" '
                f'rx="3" style="fill:var({token})">'
                f'<title>{html.escape(sample)} — {html.escape(label)}: {value:.1f}%</title>'
                f"</rect>"
            )
            if draw_w > 42:
                parts.append(
                    f'<text class="seg-label" x="{x + draw_w / 2:.1f}" y="{y + row_h / 2 + 4}" '
                    f'text-anchor="middle">{value:.0f}%</text>'
                )
            x += w

    legend = "".join(
        f'<span class="key"><span class="swatch" style="background:var({t})"></span>'
        f"{html.escape(l)}</span>"
        for _, l, t in STAGES
    )
    return f"""<section>
  <h2>Read retention</h2>
  <p class="lede">Where each sample's reads went. Segments are percentages of the
     reads entering host depletion; the final segment is what survived to
     classification.</p>
  <div class="legend">{legend}</div>
  <svg viewBox="0 0 {label_w + bar_w + 20} {height}" role="img"
       aria-label="Stacked bar chart of read retention per sample">
    {''.join(parts)}
  </svg>
</section>"""


def table(rows, columns, title, lede, empty_when=None):
    if empty_when and all(r.get(empty_when) in NA for r in rows):
        return ""
    head = "".join(f"<th>{html.escape(h)}</th>" for _, h, _ in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{fmt(row.get(key))}</td>" for key, _, fmt in columns)
        body.append(f"<tr>{cells}</tr>")
    return f"""<section>
  <h2>{html.escape(title)}</h2>
  <p class="lede">{lede}</p>
  <div class="scroll"><table><thead><tr>{head}</tr></thead>
    <tbody>{''.join(body)}</tbody></table></div>
</section>"""


def provenance_banner(record):
    if not record:
        return ('<div class="banner warn"><strong>No provenance record.</strong> '
                'This report cannot say which code produced these numbers.</div>')
    code = record.get("code", {})
    status = record.get("status", "?")
    bits = [
        f'<span><strong>Status</strong> {html.escape(status)}</span>',
        f'<span><strong>Commit</strong> <code>{html.escape(str(code.get("commit", "—"))[:12])}</code></span>',
        f'<span><strong>Nextflow</strong> {html.escape(str(record.get("environment", {}).get("nextflow_version", "—")))}</span>',
        f'<span><strong>Started</strong> {html.escape(str(record.get("run", {}).get("started", "—")))}</span>',
    ]
    meta = f'<div class="provenance">{"".join(bits)}</div>'
    if code.get("dirty"):
        n = len(code.get("modified_files", []))
        meta += (f'<div class="banner warn"><strong>Uncommitted changes.</strong> '
                 f'The working tree had {n} modified path(s), so these results are '
                 f'<em>not</em> attributable to that commit.</div>')
    return meta


# --------------------------------------------------------------------------- page

CSS = """
:root{
  --surface:#fcfcfb; --plane:#f9f9f7; --ink:#1a1a19; --ink-2:#4a4a47; --ink-3:#6f6f6a;
  --rule:#e4e4df; --accent:#2a78d6;
  --seg-1:#cde2fb; --seg-2:#3987e5; --seg-3:#0d366b; --seg-4:#1baf7a;
  --warn-bg:#fff6e5; --warn-ink:#7a4b00; --warn-rule:#fab219;
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --surface:#1a1a19; --plane:#0d0d0d; --ink:#f2f2ef; --ink-2:#c3c3bd; --ink-3:#8f8f88;
    --rule:#33332f; --accent:#3987e5;
    --seg-1:#cde2fb; --seg-2:#6da7ec; --seg-3:#256abf; --seg-4:#199e70;
    --warn-bg:#2b2410; --warn-ink:#fab219; --warn-rule:#fab219;
  }
}
:root[data-theme="dark"]{
  --surface:#1a1a19; --plane:#0d0d0d; --ink:#f2f2ef; --ink-2:#c3c3bd; --ink-3:#8f8f88;
  --rule:#33332f; --accent:#3987e5;
  --seg-1:#cde2fb; --seg-2:#6da7ec; --seg-3:#256abf; --seg-4:#199e70;
  --warn-bg:#2b2410; --warn-ink:#fab219; --warn-rule:#fab219;
}
*{box-sizing:border-box}
body{margin:0;background:var(--plane);color:var(--ink);
  font:15px/1.55 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;}
.wrap{max-width:1060px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:25px;margin:0 0 4px;letter-spacing:-.01em}
h2{font-size:17px;margin:0 0 6px;letter-spacing:-.005em}
.sub{color:var(--ink-3);margin:0 0 20px;font-size:13px}
section{background:var(--surface);border:1px solid var(--rule);border-radius:10px;
  padding:20px;margin:0 0 18px}
.lede{color:var(--ink-2);margin:0 0 14px;font-size:13.5px;max-width:74ch}
.empty{color:var(--ink-3);font-size:13.5px;margin:0}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:1px;
  background:var(--rule);padding:1px}
.tile{background:var(--surface);padding:16px 18px}
.tile-label{color:var(--ink-3);font-size:11.5px;text-transform:uppercase;
  letter-spacing:.06em}
.tile-value{font-size:27px;font-weight:600;letter-spacing:-.02em;margin:4px 0 2px;
  font-variant-numeric:tabular-nums}
.tile-note{color:var(--ink-3);font-size:12px}
.provenance{display:flex;flex-wrap:wrap;gap:18px;font-size:13px;color:var(--ink-2);
  margin:0 0 14px}
.provenance strong{color:var(--ink-3);font-weight:500;margin-right:5px}
.banner{border-radius:8px;padding:11px 14px;font-size:13.5px;margin:0 0 16px}
.banner.warn{background:var(--warn-bg);color:var(--warn-ink);
  border:1px solid var(--warn-rule)}
.legend{display:flex;flex-wrap:wrap;gap:16px;margin:0 0 14px;font-size:12.5px;
  color:var(--ink-2)}
.key{display:inline-flex;align-items:center;gap:7px}
.swatch{width:11px;height:11px;border-radius:3px;display:inline-block}
svg{width:100%;height:auto;display:block}
.cat{font-size:11.5px;fill:var(--ink-2);font-family:ui-monospace,SFMono-Regular,monospace}
.seg-label{font-size:10.5px;fill:var(--surface);font-weight:600;
  font-variant-numeric:tabular-nums;pointer-events:none}
.seg{transition:opacity .12s}
.seg:hover{opacity:.82}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:12.5px;
  font-variant-numeric:tabular-nums}
th,td{text-align:right;padding:7px 10px;border-bottom:1px solid var(--rule);
  white-space:nowrap}
th{color:var(--ink-3);font-weight:500;font-size:11.5px;text-transform:uppercase;
  letter-spacing:.05em;text-align:right}
th:first-child,td:first-child{text-align:left;
  font-family:ui-monospace,SFMono-Regular,monospace}
tbody tr:hover{background:var(--plane)}
footer{color:var(--ink-3);font-size:12px;margin-top:26px;line-height:1.7}
code{font-family:ui-monospace,SFMono-Regular,monospace;font-size:.92em}
@media (max-width:640px){.wrap{padding:20px 16px 48px}.tile-value{font-size:23px}}
"""


def validation_section(path):
    """Alignment-validation summary from validation_mask_report.tsv (absent unless it ran)."""
    if not path or not path.is_file():
        return ""
    with path.open() as fh:
        rep = list(csv.DictReader(fh, delimiter="\t"))
    if not rep:
        return ""
    tot = {k: sum(int(float(r.get(k) or 0)) for r in rep)
           for k in ("n_validated", "n_failed", "n_unsupported", "n_insufficient", "n_not_in_reference",
                     "n_unscored")}
    removed = sum(float(r.get("reads_removed") or 0) for r in rep)
    for r in rep:
        r["species"] = r["clade_name"].rsplit("|s__", 1)[-1]
    for r in rep:
        r["n_zeroed"] = int(float(r.get("n_failed") or 0)) + int(float(r.get("n_unsupported") or 0))
    worst = sorted((r for r in rep if r["n_zeroed"] > 0),
                   key=lambda r: (-r["n_zeroed"], -float(r["reads_removed"] or 0)))[:20]
    tiles = (f"<p class=\"lede\">Sample &times; species cells: <b>{tot['n_validated']:,}</b> validated, "
             f"<b>{tot['n_failed']:,}</b> failed and <b>{tot['n_unsupported']:,}</b> unsupported "
             f"(Bracken count set to 0), <b>{tot['n_insufficient']:,}</b> "
             f"insufficient (kept), <b>{tot['n_not_in_reference']:,}</b> not in the reference (kept), "
             f"<b>{tot['n_unscored']:,}</b> unscored (kept, not validated: the species was not on the validation species list). "
             f"Bracken reads removed: <b>{removed:,.0f}</b>.</p>")
    body = table(
        worst,
        [("species", "Species", html.escape),
         ("n_failed", "Samples failed", fmt_int),
         ("n_unsupported", "Unsupported", fmt_int),
         ("n_validated", "Validated", fmt_int),
         ("n_insufficient", "Insufficient", fmt_int),
         ("reads_removed", "Reads removed", fmt_int)],
        "Alignment validation: most-failed species",
        "A species <b>fails</b> in a sample when its reads are there but don't come from its genome "
        "(clumped on a few regions, low identity, or all shared with another species). It is "
        "<b>unsupported</b> when Bracken gave it at least 25 reads but fewer than 10 of them align to it. "
        "Only failed and unsupported cells are masked. Why is in "
        "<code>06_VALIDATION/validation_stats.all.tsv</code>.",
    ) if worst else "<section><h2>Alignment validation</h2><p class=\"lede\">No species failed validation in any sample.</p></section>"
    return body.replace("</h2>", "</h2>\n  " + tiles, 1)


def consensus_banner(path):
    """Say plainly which consensus mode ran; pass-through means no MetaPhlAn cross-check (A09)."""
    if not path or not path.is_file():
        return ""
    with path.open() as fh:
        st = {r["key"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}
    mode = st.get("mode", "?")
    if mode == "consensus":
        return (f'<section><h2>Consensus</h2><p class="lede">Consensus computed from '
                f'{html.escape(st.get("samples_used_for_consensus", "?"))} sample(s) with a MetaPhlAn profile; '
                f'{html.escape(st.get("genus_rows_out", "?"))} genera and {html.escape(st.get("species_rows_out", "?"))} '
                f'species kept.</p></section>')
    return (f'<section><h2>Consensus: Bracken pass-through</h2><p class="lede"><b>No MetaPhlAn cross-check '
            f'was possible</b> ({html.escape(st.get("reason", ""))}), so the Bracken tables went downstream '
            f'unfiltered. Treat taxon calls accordingly.</p></section>')


def render(rows, record, qc_path, validation_report=None, consensus_status=None) -> str:
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    filtering = table(
        rows,
        [("sample", "Sample", html.escape),
         ("library_primary_records", "Library (BAM/CRAM)", fmt_int),
         ("reads_raw", "Reads into QC", fmt_int),
         ("reads_after_fastp", "After fastp", fmt_int),
         ("pct_pass_fastp", "Pass %", fmt_pct),
         ("duplication_rate", "Dup.", lambda v: "—" if num(v) is None else f"{num(v):.3f}"),
         ("adapter_scan_applied", "Adapter scan", lambda v: html.escape(str(v or "—"))),
         ("reads_host_depleted", "Host depleted", fmt_int),
         ("pct_retained_overall", "Retained %", fmt_pct)],
        "Filtering and depletion",
        "Per-sample read accounting. <code>Adapter scan</code> records whether the "
        "234-sequence FASTA was used for that sample — <code>--adapter_trim auto</code> "
        "skips it for reads extracted from an alignment.",
    )

    classification = table(
        rows,
        [("sample", "Sample", html.escape),
         ("bracken_microbial_reads", "Microbial reads", fmt_int),
         ("cleared_consensus_gate", "Cleared gate", lambda v: html.escape(str(v or "—"))),
         ("bracken_genera", "Bracken genera", fmt_int),
         ("metaphlan_genera", "MetaPhlAn genera", fmt_int),
         ("shared_genera", "Shared genera", fmt_int),
         ("bracken_species", "Bracken species", fmt_int),
         ("metaphlan_species", "MetaPhlAn species", fmt_int),
         ("shared_species", "Shared species", fmt_int)],
        "Classification and profiler agreement",
        "<code>Microbial reads</code> is the KrakenUniq <code>root</code> clade count, i.e. the "
        "reads classified at all — the same quantity <code>--consensus_min_reads</code> gates on. A sample below "
        "the gate has no MetaPhlAn profile, so its MetaPhlAn and shared columns read "
        "&mdash;. <code>Shared</code> is the per-sample intersection of the two profilers.",
        empty_when="bracken_microbial_reads",
    )

    decontam = table(
        rows,
        [("sample", "Sample", html.escape),
         ("decontam_sample_type", "Sample type", lambda v: html.escape(str(v or "—"))),
         ("decontam_role", "Role", lambda v: html.escape(str(v or "—"))),
         ("decontam_final_status", "Status", lambda v: html.escape(str(v or "—")))],
        "Decontamination roles",
        "Which samples were treated as positives and which as negative controls. "
        "Empty unless decontamination ran.",
        empty_when="decontam_role",
    )

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>CMPipeline cohort report</title>
<style>{CSS}</style></head>
<body><div class="wrap">
<h1>CMPipeline cohort report</h1>
<p class="sub">Generated {html.escape(generated)} · source <code>{html.escape(str(qc_path.name))}</code></p>
{provenance_banner(record)}
{stat_tiles(rows)}
{retention_chart(rows)}
{filtering}
{classification}
{consensus_banner(consensus_status)}
{decontam}
{validation_section(validation_report)}
<footer>
  Every number here comes from <code>{html.escape(str(qc_path.name))}</code>; this page
  does not re-read the results tree, so the two cannot disagree.<br>
  Counts and rates only — no sequence and no clinical value. The page does carry
  per-sample identifiers, so treat it as a sample-level table.<br>
  Colour encodes ordered depletion stages; every segment is labelled and repeated in the
  tables, so nothing depends on colour alone.
</footer>
</div></body></html>"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc-summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--provenance", type=Path, default=None)
    parser.add_argument("--consensus-status", type=Path, default=None,
                        help="CONSENSUS_TAXA/consensus_status.tsv, when consensus ran")
    parser.add_argument("--validation-report", type=Path, default=None,
                        help="06_VALIDATION/validation_mask_report.tsv, when alignment validation ran")
    args = parser.parse_args()

    with args.qc_summary.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    record = None
    if args.provenance and args.provenance.is_dir():
        candidates = sorted(args.provenance.glob("*.json"),
                            key=lambda p: p.stat().st_mtime, reverse=True)
        if candidates:
            record = json.loads(candidates[0].read_text())
    elif args.provenance and args.provenance.is_file():
        record = json.loads(args.provenance.read_text())

    args.output.write_text(render(rows, record, args.qc_summary, args.validation_report, args.consensus_status))
    print(f"Cohort report: {len(rows)} sample(s) -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
