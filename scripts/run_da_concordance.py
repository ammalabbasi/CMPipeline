#!/usr/bin/env python3
"""Cross-method differential abundance concordance.

Loads results from MaAsLin, ANCOMBC, BIRDMAn, and/or Qadabra and produces a
merged concordance table with per-tool LFC/significance columns plus consensus
flags. Concordance is measured by sign agreement and rank/Pearson correlation
rather than raw LFC equality because each tool's effect size is on a different
scale.

Usage:
    python scripts/run_da_concordance.py \
        --maaslin  RESULTS/07_MAASLIN/prefix_genus_type/prefix_maaslin_full.tsv \
        --ancombc  RESULTS/06_ANCOMBC/prefix_genus_type/prefix_ancombc2_full.tsv \
        --birdman  RESULTS/08_BIRDMAN/prefix_genus_type/prefix_birdman_results.tsv \
        --qadabra  RESULTS/09_QADABRA/prefix_genus_type/concatenated_differentials.tsv \
        --level genus \
        --output   RESULTS/10_DA_CONCORDANCE/da_concordance_genus.tsv
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats


# ---------- Taxon-name extractors ---------------------------------------------

def _extract_name(taxpath, level):
    """Extract genus or species short name from a GTDB-style taxonomy path."""
    if not isinstance(taxpath, str):
        return None
    prefix = "g__" if level == "genus" else "s__"
    # Handle pipe-delimited paths
    parts = re.split(r"\|", taxpath)
    for p in reversed(parts):
        if p.startswith(prefix) and len(p) > 3:
            return p[3:]
    # Handle dot-delimited paths
    parts = re.split(r"(?<=[a-z])\.(?=[a-z]__)", taxpath)
    for p in reversed(parts):
        if p.startswith(prefix) and len(p) > 3:
            return p[3:]
    # Fallback: return as-is (may already be a short name)
    return taxpath if len(taxpath) > 0 else None


# ---------- Per-tool loaders -------------------------------------------------

def load_maaslin(path, level):
    """Load MaAsLin3/2 full results. Normalizes to (taxon, lfc, q, significant)."""
    df = pd.read_csv(path, sep="\t")

    # Normalize q-value column name across MaAsLin versions
    if "qval" not in df.columns and "qval_individual" in df.columns:
        df["qval"] = df["qval_individual"]
    if "pval" not in df.columns and "pval_individual" in df.columns:
        df["pval"] = df["pval_individual"]

    # Filter to abundance model (MaAsLin3 compound testing)
    if "model" in df.columns:
        df = df[df["model"] == "abundance"].copy()

    # Filter to Type = Tumor rows
    type_rows = df[
        df["metadata"].str.contains("Type", case=False, na=False) &
        df["value"].isin(["Tumor", "TypeTumor", "1"])
    ].copy()

    if len(type_rows) == 0:
        print(f"  WARN MaAsLin: no 'Type/Tumor' rows in {path}", flush=True)
        return None

    result = pd.DataFrame({
        "taxon_path": type_rows["feature"].values,
        "lfc": pd.to_numeric(type_rows["coef"], errors="coerce").values,
        "q":   pd.to_numeric(type_rows["qval"],  errors="coerce").values,
        "significant": (pd.to_numeric(type_rows["qval"], errors="coerce") < 0.25).values,
    })
    result["taxon"] = result["taxon_path"].apply(lambda t: _extract_name(t, level))
    result = result.dropna(subset=["taxon", "lfc"])
    result = result.groupby("taxon")[["lfc", "q", "significant"]].mean()
    result["significant"] = result["significant"] > 0.5
    return result


def load_ancombc(path, level):
    """Load ANCOMBC2 full results. Uses TypeTumor LFC and diff columns."""
    df = pd.read_csv(path, sep="\t")

    lfc_cols  = [c for c in df.columns if c.startswith("lfc_")  and "TypeTumor" in c]
    q_cols    = [c for c in df.columns if c.startswith("q_")    and "TypeTumor" in c]
    diff_cols = [c for c in df.columns if c.startswith("diff_") and "TypeTumor" in c]

    if not lfc_cols:
        print(f"  WARN ANCOMBC: no TypeTumor LFC column in {path}", flush=True)
        return None

    result = pd.DataFrame({
        "taxon_path": df["taxon"],
        "lfc": pd.to_numeric(df[lfc_cols[0]], errors="coerce"),
        "q":   pd.to_numeric(df[q_cols[0]],   errors="coerce") if q_cols else np.nan,
        "significant": (
            df[diff_cols[0]].astype(bool) if diff_cols
            else (pd.to_numeric(df[q_cols[0]], errors="coerce") < 0.05)
        ),
    })
    result["taxon"] = result["taxon_path"].apply(lambda t: _extract_name(t, level))
    result = result.dropna(subset=["taxon", "lfc"])
    result = result.groupby("taxon")[["lfc", "q", "significant"]].mean()
    result["significant"] = result["significant"] > 0.5
    return result


def load_birdman(path, level):
    """Load BIRDMAn NB results. Significance = 95% HDI excludes 0."""
    df = pd.read_csv(path, sep="\t")

    # Support both output formats
    if "lfc" not in df.columns and "coef" in df.columns:
        df = df.rename(columns={"coef": "lfc"})
    if "lfc" not in df.columns:
        print(f"  WARN BIRDMAn: no lfc column in {path}", flush=True)
        return None

    if "significant" not in df.columns and "hdi_low" in df.columns:
        df["significant"] = (df["hdi_low"] > 0) | (df["hdi_high"] < 0)

    taxon_col = next((c for c in ["taxon", "featureid", "feature", "otu"] if c in df.columns), None)
    if taxon_col is None:
        print(f"  WARN BIRDMAn: no taxon column in {path}", flush=True)
        return None

    result = pd.DataFrame({
        "taxon_path": df[taxon_col],
        "lfc": pd.to_numeric(df["lfc"], errors="coerce"),
        "q":   np.nan,
        "significant": df["significant"].astype(bool),
    })
    result["taxon"] = result["taxon_path"].apply(lambda t: _extract_name(t, level))
    result = result.dropna(subset=["taxon", "lfc"])
    result = result.groupby("taxon")[["lfc", "significant"]].mean()
    result["q"] = np.nan
    result["significant"] = result["significant"] > 0.5
    return result[["lfc", "q", "significant"]]


def load_qadabra(path, level):
    """Load Qadabra concatenated_differentials.tsv.

    The concatenated file has a 'method' column and an effect column
    (differential / ranking / log2FoldChange / lfc / coef / effect_size / logFC).
    We take the mean LFC across methods (they use different scales, so this is
    a simple aggregate); significance is majority-vote across methods.
    """
    df = pd.read_csv(path, sep="\t")

    # Find effect column
    effect_col = next(
        (c for c in ["differential", "ranking", "log2FoldChange", "lfc",
                      "coef", "effect_size", "logFC"] if c in df.columns),
        None
    )
    if effect_col is None:
        print(f"  WARN Qadabra: no effect column in {path}", flush=True)
        return None

    # Find taxon column
    taxon_col = next(
        (c for c in ["featureid", "taxon", "feature", "otu"] if c in df.columns),
        None
    )
    if taxon_col is None:
        print(f"  WARN Qadabra: no taxon column in {path}", flush=True)
        return None

    # Find significance column (optional)
    sig_col = next(
        (c for c in ["significant", "padj", "pvalue", "q_value", "qval"] if c in df.columns),
        None
    )

    df["_taxon"] = df[taxon_col].apply(lambda t: _extract_name(str(t), level))
    df["_lfc"] = pd.to_numeric(df[effect_col], errors="coerce")
    df["_sig"] = (
        df[sig_col].astype(float) < 0.05 if sig_col else np.nan
    )

    df = df.dropna(subset=["_taxon", "_lfc"])
    agg = df.groupby("_taxon").agg(
        lfc=("_lfc", "mean"),
        significant=("_sig", lambda x: x.mean() > 0.5 if not x.isna().all() else False),
    )
    agg["q"] = np.nan
    return agg[["lfc", "q", "significant"]]


LOADERS = {
    "maaslin": load_maaslin,
    "ancombc": load_ancombc,
    "birdman":  load_birdman,
    "qadabra":  load_qadabra,
}


# ---------- Main concordance logic -------------------------------------------

def build_concordance(tool_results, level):
    """Merge per-tool results into a single concordance DataFrame."""
    available = {k: v for k, v in tool_results.items() if v is not None}
    if len(available) == 0:
        raise ValueError("No tool results loaded")

    all_taxa = sorted(set.union(*[set(df.index) for df in available.values()]))
    summary = pd.DataFrame(index=all_taxa)

    for tool, df in available.items():
        tool_upper = tool.upper() if tool != "maaslin" else "MaAsLin"
        lfc_col = f"{tool_upper}_LFC"
        sig_col = f"{tool_upper}_sig"
        q_col   = f"{tool_upper}_q"
        summary[lfc_col] = df.reindex(all_taxa)["lfc"]
        summary[sig_col] = df.reindex(all_taxa)["significant"].fillna(False)
        summary[q_col]   = df.reindex(all_taxa)["q"]

    # Consensus flags
    sig_cols = [c for c in summary.columns if c.endswith("_sig")]
    lfc_cols = [c for c in summary.columns if c.endswith("_LFC")]

    summary["n_methods_sig"] = summary[sig_cols].sum(axis=1)

    signs = np.sign(summary[lfc_cols].fillna(0))
    summary["all_same_direction"] = signs.apply(
        lambda row: len(set(row[row != 0])) <= 1, axis=1
    )
    mean_lfc = summary[lfc_cols].mean(axis=1)
    summary["direction"] = np.where(mean_lfc > 0, "Enriched", "Depleted")

    n_tools = len(available)
    summary[f"consensus_{n_tools}way"] = (
        (summary["n_methods_sig"] == n_tools) & summary["all_same_direction"]
    )
    if n_tools > 2:
        summary["consensus_2way"] = (
            (summary["n_methods_sig"] >= 2) & summary["all_same_direction"]
        )

    return summary, available


def pairwise_concordance(available, lfc_cols_map):
    """Compute pairwise Pearson-r and sign-agreement among available tools."""
    tools = list(available.keys())
    rows = []
    for i, t1 in enumerate(tools):
        for t2 in tools[i+1:]:
            d1 = available[t1]["lfc"]
            d2 = available[t2]["lfc"]
            common = d1.index.intersection(d2.index)
            merged = pd.DataFrame({"a": d1.loc[common], "b": d2.loc[common]}).dropna()
            if len(merged) < 3:
                r, sign_agree = np.nan, np.nan
            else:
                r, _ = stats.pearsonr(merged["a"], merged["b"])
                sign_agree = (np.sign(merged["a"]) == np.sign(merged["b"])).mean()
            rows.append({
                "tool1": t1, "tool2": t2,
                "n_common": len(merged),
                "pearson_r": round(r, 3) if not np.isnan(r) else np.nan,
                "sign_agreement": round(sign_agree, 3) if not np.isnan(sign_agree) else np.nan,
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="DA cross-method concordance")
    parser.add_argument("--maaslin",  help="MaAsLin full results TSV")
    parser.add_argument("--ancombc",  help="ANCOMBC2 full results TSV")
    parser.add_argument("--birdman",  help="BIRDMAn results TSV")
    parser.add_argument("--qadabra",  help="Qadabra concatenated_differentials.tsv")
    parser.add_argument("--level",    default="genus", choices=["genus", "species"])
    parser.add_argument("--output",   required=True, help="Output TSV path")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    tool_paths = {
        "maaslin":  args.maaslin,
        "ancombc":  args.ancombc,
        "birdman":  args.birdman,
        "qadabra":  args.qadabra,
    }

    tool_results = {}
    for tool, path in tool_paths.items():
        if path and os.path.exists(path):
            print(f"Loading {tool}: {path}", flush=True)
            try:
                tool_results[tool] = LOADERS[tool](path, args.level)
                n = len(tool_results[tool]) if tool_results[tool] is not None else 0
                print(f"  Loaded {n} taxa", flush=True)
            except Exception as e:
                print(f"  WARN: Failed to load {tool}: {e}", flush=True)
                tool_results[tool] = None
        elif path:
            print(f"WARN: {tool} path not found: {path}", flush=True)

    tool_results = {k: v for k, v in tool_results.items() if v is not None}
    if not tool_results:
        print("ERROR: No tool results available", flush=True)
        sys.exit(1)

    print(f"\nBuilding concordance table from {list(tool_results.keys())}...", flush=True)
    summary, available = build_concordance(tool_results, args.level)
    summary.index.name = "taxon"
    summary.to_csv(args.output, sep="\t")
    print(f"Saved concordance table: {args.output}  ({len(summary)} taxa)", flush=True)

    # Pairwise stats
    pw_path = args.output.replace(".tsv", "_pairwise.tsv")
    pw = pairwise_concordance(available, {})
    pw.to_csv(pw_path, sep="\t", index=False)
    print(f"Saved pairwise concordance: {pw_path}", flush=True)

    # Print summary
    sig_cols = [c for c in summary.columns if c.endswith("_sig")]
    for col in sig_cols:
        print(f"  {col}: {summary[col].sum()} significant", flush=True)
    n_tools = len(tool_results)
    print(f"  consensus_{n_tools}way: {summary[f'consensus_{n_tools}way'].sum()}", flush=True)


if __name__ == "__main__":
    main()
