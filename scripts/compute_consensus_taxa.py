import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import sys

# Argument parser
parser = argparse.ArgumentParser(
    description="Compute consensus microbial taxa from MetaPhlAn and Bracken outputs, "
                "filter common taxa based on samples with at least 100,000 reads, and generate a visualization."
)

parser.add_argument("--metaphlan", required=True, help="Path to the MetaPhlAn merged abundance table (genus level).")
parser.add_argument("--bracken_genus", required=True, help="Path to the Bracken genus MPA report.")
parser.add_argument("--bracken_species", required=True, help="Path to the Bracken species MPA report.")
parser.add_argument("--output", required=True, help="Path to save the output visualization (PDF).")
parser.add_argument("--output_common_genus",required=True, help="Path to save filtered common genus report (TXT).")
parser.add_argument("--output_common_species", required=True, help="Path to save filtered common species report (TXT).")

args = parser.parse_args()


# Load MetaPhlAn data
metaphlan_df = pd.read_csv(args.metaphlan, sep="\t", header=1)
# Extract genus name from MetaPhlAn clade_name (format: g__GenusName)
metaphlan_df['genus_name'] = metaphlan_df['clade_name'].str.replace(r'^g__', '', regex=True)
metaphlan_df = metaphlan_df.loc[~metaphlan_df['genus_name'].str.contains('GGB|_unclassified', na=False)]
metaphlan_df = metaphlan_df.dropna(subset=['genus_name'])
metaphlan_sample_cols = [c for c in metaphlan_df.columns if c not in ('clade_name', 'genus_name')]
print(f"MetaPhlAn: {len(metaphlan_df)} genera, {len(metaphlan_sample_cols)} samples")


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
bracken_genus_df = bracken_genus_df.loc[bracken_genus_df['genus_name'] != 'Cutibacterium']
print(f"Bracken genus: {len(bracken_genus_df)} genera, {len(bracken_sample_cols)} samples")

# Filter to samples above read count threshold (using Bracken counts)
sample_sums = bracken_genus_df[bracken_sample_cols].sum()
samples_above_threshold = sample_sums[sample_sums >= 100000].index.tolist()
print(f"Samples above 100K threshold: {len(samples_above_threshold)}")

if len(samples_above_threshold) < 2:
    print("Consensus not possible: fewer than 2 samples have at least 100,000 reads.")
    sys.exit(0)


# Load Bracken species data
bracken_species_df = pd.read_csv(args.bracken_species, sep="\t")
bracken_species_df.rename(columns={'#Classification': 'clade_name'}, inplace=True)
bracken_species_df = bracken_species_df[~bracken_species_df['clade_name'].str.contains('Viruses', na=False)]
bracken_species_df.columns = bracken_species_df.columns.str.replace(r'\.bracken\.S\.krakenreport\.txt$', '', regex=True)
bracken_species_df = bracken_species_df.dropna(subset=['clade_name'])
bracken_species_df['genus_name'] = bracken_species_df['clade_name'].str.split(r'\|g__', expand=True)[1].str.split(r'\|s__', expand=True)[0]
bracken_species_df['species_name'] = bracken_species_df['clade_name'].str.split(r'\|g__', expand=True)[1].str.split(r'\|s__', expand=True)[1]
bracken_species_df = bracken_species_df.dropna(subset=['species_name'])
bracken_species_df = bracken_species_df.loc[bracken_species_df['genus_name'] != 'Cutibacterium']


# Calculate proportions independently for each tool using genus_name as key
# Find common samples between MetaPhlAn and Bracken for fair comparison
common_samples = sorted(set(samples_above_threshold) & set(metaphlan_sample_cols))
print(f"Common samples for consensus: {len(common_samples)}")

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
common_genera = common_taxa['genus_name'].tolist()
bracken_genus_df = bracken_genus_df[bracken_genus_df['genus_name'].isin(common_genera)]
bracken_species_df = bracken_species_df[bracken_species_df['genus_name'].isin(common_genera)]
print(f"Bracken genus rows after filter: {len(bracken_genus_df)}")
print(f"Bracken species rows after filter: {len(bracken_species_df)}")

# Save filtered Bracken genus/species data
bracken_genus_df.to_csv(args.output_common_genus, sep="\t", index=False)
bracken_species_df.to_csv(args.output_common_species, sep="\t", index=False)

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
