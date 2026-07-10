#!/usr/bin/env python3
"""BIRDMAn Negative Binomial differential abundance analysis.

Runs the birdman NegativeBinomialSingle model per-feature over a corrected
count table. Depth (sequencing effort) is handled internally by the model as
log(sum of counts per sample) — equivalent to TSS-based depth correction in
linear models.

Usage:
    conda run -n q2-birdman python scripts/run_birdman.py \
        --otu RESULTS/05_BATCH_CORRECTION/corrected/ConQuR_tuned.tsv \
        --meta metadata.tsv \
        --prefix results/birdman_genus \
        --formula "Type" \
        --level genus
"""

import argparse
import os
import sys
import warnings

import arviz as az
import numpy as np
import pandas as pd
from biom.table import Table as BiomTable

warnings.filterwarnings("ignore")


def filter_to_level(df, level):
    if level == "genus":
        cols = [c for c in df.columns if c.split("|")[-1].startswith("g__")]
    elif level == "species":
        cols = [c for c in df.columns if "|s__" in c]
    else:
        raise ValueError(f"Unknown level: {level}")
    return df[cols]


def extract_short_name(taxpath, level):
    parts = taxpath.split("|")
    prefix = "g__" if level == "genus" else "s__"
    for p in reversed(parts):
        if p.startswith(prefix) and len(p) > 3:
            return p[3:]
    return taxpath


def parse_args():
    parser = argparse.ArgumentParser(description="BIRDMAn NB differential abundance")
    parser.add_argument("--otu", required=True, help="Corrected count TSV (samples x taxa)")
    parser.add_argument("--meta", required=True, help="Metadata TSV")
    parser.add_argument("--prefix", required=True, help="Output prefix")
    parser.add_argument("--formula", default="Type", help="Patsy formula (default: Type)")
    parser.add_argument("--level", default="genus", choices=["genus", "species"])
    parser.add_argument("--reference", default="Normal",
                        help="Reference level of the Type variable (default: Normal)")
    parser.add_argument("--target", default="Tumor",
                        help="Target level of the Type variable (default: Tumor)")
    parser.add_argument("--min_prevalence", type=float, default=0.10,
                        help="Min prevalence fraction to keep a taxon (default: 0.10)")
    parser.add_argument("--num_iter", type=int, default=500)
    parser.add_argument("--num_warmup", type=int, default=None)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def load_table(otu_path, meta_path, level, min_prevalence):
    df = pd.read_csv(otu_path, sep="\t", index_col=0)
    df.index = df.index.astype(str)

    df = filter_to_level(df, level)
    print(f"  After level filter: {df.shape[1]} taxa", flush=True)

    nonzero = df.sum(axis=0) > 0
    df = df.loc[:, nonzero]

    prev = (df > 0).mean(axis=0)
    df = df.loc[:, prev >= min_prevalence]
    print(f"  After prevalence filter (>={min_prevalence}): {df.shape[1]} taxa", flush=True)

    # Load metadata
    meta = pd.read_csv(meta_path, sep="\t")
    id_candidates = ["donor_id", "sample_id", "sampleid", "Sample_ID",
                     "SampleID", "patient", "Patient", "#SampleID", "sample-id"]
    id_col = next((c for c in id_candidates if c in meta.columns), None)
    if id_col is None:
        raise ValueError("Cannot find sample ID column in metadata")
    meta = meta.set_index(id_col)
    meta.index = meta.index.astype(str)

    common = sorted(set(df.index) & set(meta.index))
    print(f"  Common samples: {len(common)}", flush=True)
    df = df.loc[common]
    meta = meta.loc[common].copy()

    # Round to integer counts (ConQuR output is float on count scale)
    count_matrix = df.T.values.astype(int)  # taxa x samples for BIOM

    table = BiomTable(
        count_matrix,
        observation_ids=df.columns.tolist(),
        sample_ids=df.index.tolist(),
    )
    return table, meta


def fit_birdman(table, meta, formula, args):
    import birdman

    coef_name = f"C(Type, Treatment('{args.reference}'))[T.{args.target}]"

    results = []
    feature_ids = list(table.ids(axis="observation"))
    n = len(feature_ids)
    print(f"  Fitting {n} features with NegativeBinomial...", flush=True)

    for i, fid in enumerate(feature_ids):
        if (i + 1) % 50 == 0 or i == 0:
            print(f"  Feature {i+1}/{n}: {fid[:60]}", flush=True)
        try:
            model = birdman.NegativeBinomialSingle(
                table=table,
                feature_id=fid,
                formula=formula,
                metadata=meta,
                num_iter=args.num_iter,
                num_warmup=args.num_warmup,
                chains=args.chains,
                seed=args.seed,
            )
            model.compile_model()
            model.fit_model()
            inf = model.to_inference()

            summary = az.summary(
                inf,
                var_names=["beta_var"],
                hdi_prob=0.95,
            )

            if coef_name in summary.index:
                row = summary.loc[coef_name]
                lfc = float(row["mean"])
                hdi_lo = float(row["hdi_2.5%"])
                hdi_hi = float(row["hdi_97.5%"])
                significant = (hdi_lo > 0) or (hdi_hi < 0)
            else:
                # Formula covariate name may differ; take the first non-intercept covariate
                non_intercept = [idx for idx in summary.index if "Intercept" not in idx]
                if non_intercept:
                    row = summary.loc[non_intercept[0]]
                    lfc = float(row["mean"])
                    hdi_lo = float(row["hdi_2.5%"])
                    hdi_hi = float(row["hdi_97.5%"])
                    significant = (hdi_lo > 0) or (hdi_hi < 0)
                    coef_name = non_intercept[0]
                else:
                    continue

            results.append({
                "taxon": fid,
                "taxon_short": extract_short_name(fid, args.level),
                "coef_name": coef_name,
                "lfc": lfc,
                "hdi_low": hdi_lo,
                "hdi_high": hdi_hi,
                "significant": significant,
            })
        except Exception as e:
            print(f"  WARN: Feature {fid[:60]} failed: {e}", flush=True)

    return pd.DataFrame(results)


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.prefix) or ".", exist_ok=True)

    print(f"Loading table (level={args.level}, min_prevalence={args.min_prevalence})...", flush=True)
    table, meta = load_table(args.otu, args.meta, args.level, args.min_prevalence)

    print(f"\nFitting BIRDMAn NB (formula={args.formula})...", flush=True)
    results = fit_birdman(table, meta, args.formula, args)

    if results.empty:
        print("ERROR: No results produced", flush=True)
        sys.exit(1)

    n_sig = results["significant"].sum()
    print(f"\nResults: {len(results)} features, {n_sig} significant (95% HDI excludes 0)", flush=True)

    out_path = f"{args.prefix}_birdman_results.tsv"
    results.to_csv(out_path, sep="\t", index=False)
    print(f"Saved: {out_path}", flush=True)

    # Summary
    summary_path = f"{args.prefix}_birdman_summary.txt"
    with open(summary_path, "w") as f:
        f.write("BIRDMAn NB Analysis Summary\n")
        f.write("============================\n\n")
        f.write(f"Formula: {args.formula}\n")
        f.write(f"Level: {args.level}\n")
        f.write(f"Samples: {table.shape[1]}\n")
        f.write(f"Features fitted: {len(results)}\n")
        f.write(f"Significant (95% HDI excludes 0): {n_sig}\n")
        f.write(f"Reference level: {args.reference}\n")
        f.write(f"Target level: {args.target}\n")
    print(f"Saved: {summary_path}", flush=True)


if __name__ == "__main__":
    main()
