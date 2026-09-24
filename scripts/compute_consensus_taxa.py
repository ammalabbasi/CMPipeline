import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import sys

# Argument parser
parser = argparse.ArgumentParser(
    description="Compute consensus microbial taxa from MetaPhlAn and Bracken outputs, learning from the "
                "samples MetaPhlAn ran on (those that cleared the KrakenUniq root-count gate upstream; "
                "decision 18), and generate a visualization."
)

parser.add_argument("--metaphlan", default=None,
                    help="Path to the MetaPhlAn merged abundance table (genus level). Omit when no "
                         "sample cleared --min_reads: MetaPhlAn then never ran, and Bracken passes through.")
parser.add_argument("--bracken_genus", required=True, help="Path to the Bracken genus MPA report.")
parser.add_argument("--bracken_species", required=True, help="Path to the Bracken species MPA report.")
parser.add_argument("--min_reads", type=float, default=100000,
                    help="Reported only (samples_above_min_reads in consensus_status.tsv): the gate itself is "
                         "applied upstream on the KrakenUniq root count, and consensus uses every sample "
                         "MetaPhlAn ran on (default: 100000).")
parser.add_argument("--output", required=True, help="Path to save the output visualization (PDF).")
parser.add_argument("--output_common_genus",required=True, help="Path to save filtered common genus report (TXT).")
parser.add_argument("--output_common_species", required=True, help="Path to save filtered common species report (TXT).")
parser.add_argument("--output_library_sizes", default="bracken_library_sizes.tsv",
                    help="Path to save each sample's total Bracken reads (genus and species tables, "
                         "before any consensus filtering). Decontamination uses it to tell a sample "
                         "that is all-zero after filtering from an empty library.")
parser.add_argument("--output_status", default="consensus_status.tsv",
                    help="Path to save which path ran (consensus or bracken_passthrough) and why, as counts only.")

args = parser.parse_args()


# Load MetaPhlAn data, if MetaPhlAn ran on any sample at all
if args.metaphlan:
    metaphlan_df = pd.read_csv(args.metaphlan, sep="\t", header=1)
    # Extract genus name from MetaPhlAn clade_name (format: g__GenusName)
    metaphlan_df['genus_name'] = metaphlan_df['clade_name'].str.replace(r'^g__', '', regex=True)
    metaphlan_df = metaphlan_df.loc[~metaphlan_df['genus_name'].str.contains('GGB|_unclassified', na=False)]
    metaphlan_df = metaphlan_df.dropna(subset=['genus_name'])
    metaphlan_sample_cols = [c for c in metaphlan_df.columns if c not in ('clade_name', 'genus_name')]
    print(f"MetaPhlAn: {len(metaphlan_df)} genera, {len(metaphlan_sample_cols)} samples")
else:
    metaphlan_df = None
    metaphlan_sample_cols = []
    print("MetaPhlAn: no table (no sample cleared the read threshold)")


# Load Bracken genus data
bracken_genus_df = pd.read_csv(args.bracken_genus, sep="\t")
bracken_genus_df.rename(columns={'#Classification': 'clade_name'}, inplace=True)
bracken_genus_df = bracken_genus_df[~bracken_genus_df['clade_name'].str.contains(r'Viruses|\|s__', na=False)]
bracken_genus_df.columns = bracken_genus_df.columns.str.replace(r'\.bracken\.G\.krakenreport\.txt$', '', regex=True)
bracken_sample_cols = [c for c in bracken_genus_df.columns if c != 'clade_name']
bracken_genus_df = bracken_genus_df.dropna(subset=['clade_name'])
# Extract genus name from full taxonomy (d__...|g__GenusName)
bracken_genus_df['genus_name'] = bracken_genus_df['clade_name'].str.split(r'\|g__', expand=True)[1]
bracken_genus_df = bracken_genus_df.dropna(subset=['genus_name'])
print(f"Bracken genus: {len(bracken_genus_df)} genera, {len(bracken_sample_cols)} samples")

# Cohort-level eligibility: calculate totals across the merged Bracken genus table,
# then use the same eligible sample set for both Bracken and MetaPhlAn prevalence.
sample_sums = bracken_genus_df[bracken_sample_cols].sum()
samples_above_threshold = sample_sums[sample_sums >= args.min_reads].index.tolist()
print(f"Samples above {args.min_reads:g}-read threshold: {len(samples_above_threshold)}")


# Load Bracken species data
bracken_species_df = pd.read_csv(args.bracken_species, sep="\t")
bracken_species_df.rename(columns={'#Classification': 'clade_name'}, inplace=True)
bracken_species_df = bracken_species_df[~bracken_species_df['clade_name'].str.contains('Viruses', na=False)]
bracken_species_df.columns = bracken_species_df.columns.str.replace(r'\.bracken\.S\.krakenreport\.txt$', '', regex=True)
bracken_species_df = bracken_species_df.dropna(subset=['clade_name'])
bracken_species_df['genus_name'] = bracken_species_df['clade_name'].str.split(r'\|g__', expand=True)[1].str.split(r'\|s__', expand=True)[0]
bracken_species_df['species_name'] = bracken_species_df['clade_name'].str.split(r'\|g__', expand=True)[1].str.split(r'\|s__', expand=True)[1]
bracken_species_df = bracken_species_df.dropna(subset=['species_name'])

# Per-sample library sizes, from the cleaned Bracken tables BEFORE any consensus filtering.
species_sample_cols = [c for c in bracken_species_df.columns
                       if c not in ('clade_name', 'genus_name', 'species_name')]
genus_totals = bracken_genus_df[bracken_sample_cols].apply(pd.to_numeric, errors='coerce').sum()
species_totals = bracken_species_df[species_sample_cols].apply(pd.to_numeric, errors='coerce').sum()
pd.DataFrame({
    "sample": bracken_sample_cols,
    "bracken_genus_total": [genus_totals.get(c, 0) for c in bracken_sample_cols],
    "bracken_species_total": [species_totals.get(c, 0) for c in bracken_sample_cols],
}).to_csv(args.output_library_sizes, sep="\t", index=False)
# Calculate proportions independently for each tool using genus_name as key
# Find common samples between MetaPhlAn and Bracken for fair comparison
# Consensus learns from exactly the samples MetaPhlAn ran on. The Bracken module already gated
# them on --consensus_min_reads using the KrakenUniq ROOT count; re-thresholding here on Bracken
# genus totals (a different number) could drop a gated sample and even flip the run to
# pass-through (audit A05; decided 2026-09-23). samples_above_threshold is reported, not used.
common_samples = sorted(set(metaphlan_sample_cols) & set(bracken_sample_cols))
print(f"Common samples for consensus: {len(common_samples)}")

# Consensus is a cohort-level genus filter learned from the samples deep enough for
# MetaPhlAn. With fewer than 2 such samples there is nothing to learn it from, so the
# Bracken tables pass through unfiltered (same cleaning, same file names) and the
# status file records that this happened. Samples below --min_reads never feed the
# consensus; their Bracken counts always go forward, filtered or not.
passthrough = len(common_samples) < 2
if passthrough:
    reason = (f"{len(metaphlan_sample_cols)} sample(s) cleared the MetaPhlAn gate "
              f"(>= {args.min_reads:g} KrakenUniq root reads), "
              f"{len(common_samples)} also profiled by MetaPhlAn; consensus needs >= 2")
    print(f"WARNING: consensus skipped, Bracken passed through unfiltered: {reason}")

if not passthrough:
    # Bracken proportions: fraction of samples with non-zero counts per genus
    bracken_props = {}
    for _, row in bracken_genus_df.iterrows():
        genus = row['genus_name']
        vals = row[common_samples].values.astype(float)
        bracken_props[genus] = (vals > 0).sum() / len(common_samples)

    # MetaPhlAn proportions: fraction of samples with non-zero abundance per genus
    metaphlan_props = {}
    for _, row in metaphlan_df.iterrows():
        genus = row['genus_name']
        vals = row[common_samples].values.astype(float) if all(s in metaphlan_df.columns for s in common_samples) else []
        if len(vals) > 0:
            metaphlan_props[genus] = (vals > 0).sum() / len(common_samples)

    # Build proportions DataFrame
    all_genera = sorted(set(bracken_props.keys()) | set(metaphlan_props.keys()))
    proportions = []
    for genus in all_genera:
        proportions.append([genus, bracken_props.get(genus, 0), metaphlan_props.get(genus, 0)])

    proportions_df = pd.DataFrame(proportions, columns=['genus_name', 'Bracken', 'Metaphlan'])
    common_taxa = proportions_df[(proportions_df['Bracken'] >= 0.01) & (proportions_df['Metaphlan'] >= 0.01)]
    print(f"Total genera assessed: {len(proportions_df)}")
    print(f"Common genera (>=1% prevalence in both): {len(common_taxa)}")


# Filter all samples based on common genera (retain all samples, not just filtered ones)
if not passthrough:
    common_genera = common_taxa['genus_name'].tolist()
    bracken_genus_df = bracken_genus_df[bracken_genus_df['genus_name'].isin(common_genera)]
    bracken_species_df = bracken_species_df[bracken_species_df['genus_name'].isin(common_genera)]
print(f"Bracken genus rows after filter: {len(bracken_genus_df)}")
print(f"Bracken species rows after filter: {len(bracken_species_df)}")

# Save count-only tables. `genus_name`/`species_name` are helper columns, not samples,
# and the decontamination script treats every column but `clade_name` as a sample.
bracken_genus_output = bracken_genus_df.drop(columns=['genus_name'])
bracken_species_output = bracken_species_df.drop(columns=['genus_name', 'species_name'])
bracken_genus_output.to_csv(args.output_common_genus, sep="\t", index=False)
bracken_species_output.to_csv(args.output_common_species, sep="\t", index=False)

# Counts only -- no sample identifiers, so this is safe to quote in a report.
status = [
    ("mode", "bracken_passthrough" if passthrough else "consensus"),
    ("reason", reason if passthrough else "consensus computed"),
    ("min_reads", f"{args.min_reads:g}"),
    ("samples_in_bracken", len(bracken_sample_cols)),
    ("samples_above_min_reads", len(samples_above_threshold)),
    ("samples_with_metaphlan", len(metaphlan_sample_cols)),
    ("samples_used_for_consensus", len(common_samples)),
    ("consensus_eligibility_basis", "samples with a MetaPhlAn profile (Bracken-module gate on KrakenUniq root count)"),
    ("genus_rows_out", len(bracken_genus_output)),
    ("species_rows_out", len(bracken_species_output)),
]
pd.DataFrame(status, columns=["key", "value"]).to_csv(args.output_status, sep="\t", index=False)

if passthrough:
    # Keep the declared PDF output; say plainly why there is no scatter.
    plt.figure(figsize=(6, 6))
    plt.axis('off')
    plt.text(0.5, 0.5, "Consensus skipped - Bracken passed through unfiltered\n\n" + reason,
             ha='center', va='center', wrap=True, fontsize=10)
    plt.savefig(args.output)
    sys.exit(0)

# Plot Results
plt.figure(figsize=(6, 6))
sns.scatterplot(data=proportions_df, x='Bracken', y='Metaphlan', s=60, alpha=0.7, color='lightgray', edgecolor='k', linewidth=0.5)
sns.scatterplot(data=common_taxa, x='Bracken', y='Metaphlan', s=60, color='firebrick', edgecolor='k', linewidth=0.5)

texts = [plt.text(row['Bracken'], row['Metaphlan'], row['genus_name'], color='firebrick', fontsize=7)
         for _, row in common_taxa.iterrows()]
adjust_text(texts, arrowprops=dict(arrowstyle="->", color='k', lw=1))

plt.axvline(x=0.01, color='gray', linestyle='--', linewidth=0.7)
plt.axhline(y=0.01, color='gray', linestyle='--', linewidth=0.7)
plt.xlabel('Proportion of samples\n Bracken', fontsize=12)
plt.ylabel('Proportion of samples\n MetaPhlan', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(args.output, dpi=300)
plt.show()
