# -*- coding: utf-8 -*-
"""
Combined PFAM and KEGG enrichment plots into one SVG with panels A and B
Created on: 2025-10-04
@author: hayat 
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import Table2x2
import matplotlib.pyplot as plt
import seaborn as sns


plt.rcParams["svg.fonttype"] = "none"  
plt.rcParams["font.family"] = "Arial"  
plt.rcParams["font.size"] = 10
# -----------------------------
# File paths
# -----------------------------
PLASMID_LIST_FILE = r"C:/Users/hayat/Downloads/R_files/data/plasmids_in_Agelas.tsv"  # first column: plasmid IDs
EGGNOG_ANNOTATIONS = r"C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations"
ANNOTATIONS_KEGG_CSV = r"C:/Users/hayat/Downloads/R_files/data/agelas_enriched_pathways_annotations.csv"
OUTPUT_SVG = r"C:/Users/hayat/Downloads/R_files/graphs/agelas_enrichment_panels_AB.svg"

# -----------------------------
# Helper functions
# -----------------------------

def find_header_line(path):
    with open(path, 'r', encoding='utf-8') as fh:
        for i, line in enumerate(fh):
            if line.startswith("#query"):
                return i
    raise ValueError("Could not find header line starting with '#query' in eggNOG file")


def safe_split(pathway):
    if isinstance(pathway, str):
        return [item.strip() for item in pathway.split(',')]
    return []


def unify_kegg_id(kegg_id):
    if not isinstance(kegg_id, str):
        return kegg_id
    prefixes = ['ko', 'map', 'rn', 'ec']
    for prefix in prefixes:
        if kegg_id.startswith(prefix):
            num_part = kegg_id[len(prefix):]
            num_part_padded = num_part.zfill(5)
            return 'ko' + num_part_padded
    return kegg_id

# PFAM group mapping 
pfam_group_map = {
    'Arm-DNA-bind_4': 'Phage integrases',
    'Phage_int_SAM_3': 'Phage integrases',
    'Phage_integrase': 'Phage integrases',
    'Arm-DNA-bind_3': 'Phage integrases',
    'GFO_IDH_MocA': 'GFO_IDH_MocA_C',
    'GFO_IDH_MocA_C': 'GFO_IDH_MocA_C',
}


def map_to_group(pfam):
    return pfam_group_map.get(pfam, pfam)

# -----------------------------
# Load data
# -----------------------------
print("Loading data...")
plasmid_df = pd.read_csv(PLASMID_LIST_FILE, sep='\t', header=None)
plasmids = set(plasmid_df.iloc[:, 0].astype(str))

header_line = find_header_line(EGGNOG_ANNOTATIONS)
eggnog_df = pd.read_csv(EGGNOG_ANNOTATIONS, sep='\t', skiprows=header_line, low_memory=False, dtype=str)

# Extract plasmid id from '#query' like 'DRR066339_127ctg_1' -> 'DRR066339_127ctg'
eggnog_df['Plasmid'] = eggnog_df['#query'].str.extract(r'(^[^_]+_[^_]+)')

# Ensure expected columns exist (PFAMs, KEGG_Pathway)
for col in ['PFAMs', 'KEGG_Pathway']:
    if col not in eggnog_df.columns:
        eggnog_df[col] = '-'

# Split foreground/background once
foreground = eggnog_df[eggnog_df['Plasmid'].isin(plasmids)].copy()
background = eggnog_df[~eggnog_df['Plasmid'].isin(plasmids)].copy()

# -----------------------------
# PFAM enrichment (Panel A)
# -----------------------------
print("Running PFAM enrichment...")
# Remove missing PFAMs and explode
fg_pfam = foreground[foreground['PFAMs'] != '-'].copy()
bg_pfam = background[background['PFAMs'] != '-'].copy()

fg_pfam['PFAMs_list'] = fg_pfam['PFAMs'].str.split(',')
bg_pfam['PFAMs_list'] = bg_pfam['PFAMs'].str.split(',')

fg_exp = fg_pfam.explode('PFAMs_list')
bg_exp = bg_pfam.explode('PFAMs_list')

fg_exp['PFAMs_list'] = fg_exp['PFAMs_list'].str.strip()
bg_exp['PFAMs_list'] = bg_exp['PFAMs_list'].str.strip()

fg_exp['PFAM_Group'] = fg_exp['PFAMs_list'].apply(map_to_group)
bg_exp['PFAM_Group'] = bg_exp['PFAMs_list'].apply(map_to_group)

fg_counts = fg_exp['PFAM_Group'].value_counts()
bg_counts = bg_exp['PFAM_Group'].value_counts()
all_pfams = fg_counts.index.union(bg_counts.index)

fg_total = len(fg_exp)
bg_total = len(bg_exp)

pfam_results = []
for pfam in all_pfams:
    fg_hits = int(fg_counts.get(pfam, 0))
    fg_miss = fg_total - fg_hits
    bg_hits = int(bg_counts.get(pfam, 0))
    bg_miss = bg_total - bg_hits

    contingency_table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    try:
        table = Table2x2(contingency_table)
        odds_ratio = table.oddsratio
        ci_low, ci_upp = table.oddsratio_confint()
    except Exception:
        # fallback when Table2x2 fails (e.g., zeros)
        odds_ratio, _ = fisher_exact(contingency_table, alternative='greater')
        ci_low, ci_upp = (np.nan, np.nan)

    p_value = fisher_exact(contingency_table, alternative='greater')[1]
    pfam_results.append({
        'PFAMs': pfam,
        'Foreground_Count': fg_hits,
        'Background_Count': bg_hits,
        'Odds_Ratio': odds_ratio,
        'CI_Lower': ci_low,
        'CI_Upper': ci_upp,
        'P_Value': p_value
    })

pfam_df = pd.DataFrame(pfam_results)
pfam_df['P_Value_Adjusted'] = multipletests(pfam_df['P_Value'], method='fdr_bh')[1]

# Apply thresholds (tweak if needed)
pfam_sig = pfam_df[(pfam_df['P_Value_Adjusted'] < 0.05) & (pfam_df['Odds_Ratio'] > 8) & (pfam_df['Foreground_Count'] >= 10)].copy()
if pfam_sig.empty:
    print("Warning: no PFAM groups passed the significance filters. Consider lowering thresholds.")

# Prepare plotting columns
pfam_sig = pfam_sig.replace([np.inf, -np.inf], np.nan).dropna(subset=['Odds_Ratio'])
pfam_sig['log2_OR'] = np.log2(pfam_sig['Odds_Ratio'])
pfam_sig['log2_CI_Lower'] = np.log2(pfam_sig['CI_Lower'])
pfam_sig['log2_CI_Upper'] = np.log2(pfam_sig['CI_Upper'])
pfam_sig.sort_values('log2_OR', inplace=True)

# -----------------------------
# KEGG enrichment (Panel B)
# -----------------------------
print("Running KEGG enrichment...")
fg_kegg = foreground[foreground['KEGG_Pathway'] != '-'].copy()
bg_kegg = background[background['KEGG_Pathway'] != '-'].copy()

fg_kegg['KEGG_Pathway'] = fg_kegg['KEGG_Pathway'].apply(safe_split)
bg_kegg['KEGG_Pathway'] = bg_kegg['KEGG_Pathway'].apply(safe_split)

fg_kegg = fg_kegg.explode('KEGG_Pathway')
bg_kegg = bg_kegg.explode('KEGG_Pathway')

fg_kegg['KEGG_Pathway_Unified'] = fg_kegg['KEGG_Pathway'].apply(unify_kegg_id)
bg_kegg['KEGG_Pathway_Unified'] = bg_kegg['KEGG_Pathway'].apply(unify_kegg_id)

fg_counts_k = fg_kegg['KEGG_Pathway_Unified'].value_counts()
bg_counts_k = bg_kegg['KEGG_Pathway_Unified'].value_counts()
all_pathways = fg_counts_k.index.union(bg_counts_k.index)

fg_total_k = len(fg_kegg)
bg_total_k = len(bg_kegg)

kegg_results = []
for pathway in all_pathways:
    fg_hits = int(fg_counts_k.get(pathway, 0))
    fg_miss = fg_total_k - fg_hits
    bg_hits = int(bg_counts_k.get(pathway, 0))
    bg_miss = bg_total_k - bg_hits

    contingency_table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')

    kegg_results.append({
        'KEGG_Pathway': pathway,
        'Foreground_Count': fg_hits,
        'Background_Count': bg_hits,
        'Odds_Ratio': odds_ratio,
        'P_Value': p_value
    })

kegg_df = pd.DataFrame(kegg_results)
if not kegg_df.empty:
    kegg_df['P_Value_Adjusted'] = multipletests(kegg_df['P_Value'], method='fdr_bh')[1]
else:
    kegg_df['P_Value_Adjusted'] = []

# thresholds
kegg_sig = kegg_df[(kegg_df['P_Value_Adjusted'] < 0.05) & (kegg_df['Odds_Ratio'] > 1.0) & (kegg_df['Foreground_Count'] >= 20)].copy()
if kegg_sig.empty:
    print("Warning: no KEGG pathways passed the significance filters. Consider lowering thresholds.")

# compute log2 OR and approximate 95% CI on log2 scale (add pseudo-counts for stability)

if not kegg_sig.empty:
    def compute_log2_ci(row):
        a = max(row['Foreground_Count'], 1)
        b = max(fg_total_k - row['Foreground_Count'], 1)
        c = max(row['Background_Count'], 1)
        d = max(bg_total_k - row['Background_Count'], 1)
        # avoid zero or infinite OR
        log_or = np.log2(row['Odds_Ratio']) if row['Odds_Ratio'] > 0 else 0.0
        se = np.sqrt(1/a + 1/b + 1/c + 1/d) / np.log(2)
        ci_lower = log_or - 1.96 * se
        ci_upper = log_or + 1.96 * se
        return pd.Series([log_or, ci_lower, ci_upper])

    kegg_sig[['log2_OR', 'log2_CI_Lower', 'log2_CI_Upper']] = kegg_sig.apply(compute_log2_ci, axis=1)
    # add Y labels from annotation CSV if present
    try:
        annotations = pd.read_csv(ANNOTATIONS_KEGG_CSV, header=None, names=['ID', 'Label'])
        annotations['ID'] = annotations['ID'].str.strip()
        annotations['Label'] = annotations['Label'].str.strip()
        kegg_sig = kegg_sig.merge(annotations, left_on='KEGG_Pathway', right_on='ID', how='left')
        kegg_sig['Y_Label'] = kegg_sig.apply(lambda row: f"{row['Label']} ({row['KEGG_Pathway']})" if pd.notna(row.get('Label')) else row['KEGG_Pathway'], axis=1)
    except Exception:
        kegg_sig['Y_Label'] = kegg_sig['KEGG_Pathway']
    kegg_sig.sort_values('log2_OR', inplace=True)

# -----------------------------
# Plot combined figure (Panel A: PFAM, Panel B: KEGG)
# -----------------------------
print("Creating combined figure and saving to SVG...")

sns.set(style='whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={'width_ratios': [1, 1]})

# Panel A: PFAM
ax = axes[0]
ax.set_title('(a)', loc='left', fontsize=12, fontweight='bold')
if not pfam_sig.empty:
    y_pos = np.arange(len(pfam_sig))
    ax.errorbar(
        x=pfam_sig['log2_OR'],
        y=y_pos,
        xerr=[pfam_sig['log2_OR'] - pfam_sig['log2_CI_Lower'], pfam_sig['log2_CI_Upper'] - pfam_sig['log2_OR']],
        fmt='o', ecolor='gray', elinewidth=1.5, capsize=3
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pfam_sig['PFAMs'])
else:
    ax.text(0.5, 0.5, 'No PFAMs passed filters', ha='center', va='center')
ax.axvline(0, color='red', linestyle='--')
ax.set_xlabel(r'$\log_2(\mathrm{Odds\ Ratio})$')
ax.set_ylabel('PFAM domain group')

# Panel B: KEGG
ax2 = axes[1]
ax2.set_title('(b)', loc='left', fontsize=12, fontweight='bold')
if not kegg_sig.empty:
    y_pos2 = np.arange(len(kegg_sig))
    ax2.errorbar(
        x=kegg_sig['log2_OR'],
        y=y_pos2,
        xerr=[kegg_sig['log2_OR'] - kegg_sig['log2_CI_Lower'], kegg_sig['log2_CI_Upper'] - kegg_sig['log2_OR']],
        fmt='o', ecolor='gray', elinewidth=1.5, capsize=3
    )
    ax2.set_yticks(y_pos2)
    ax2.set_yticklabels(kegg_sig['Y_Label'])
else:
    ax2.text(0.5, 0.5, 'No KEGG pathways passed filters', ha='center', va='center')
ax2.axvline(0, color='red', linestyle='--')
ax2.set_xlabel(r'$\log_2(\mathrm{Odds\ Ratio})$')
ax2.set_ylabel('')

plt.tight_layout()
fig.savefig(OUTPUT_SVG, format='svg')
print(f"Saved combined figure to: {OUTPUT_SVG}")

# Optionally also show the figure when running interactively
# plt.show()

# Save intermediate tables (uncomment to enable)
# pfam_sig.to_csv(r"C:/Users/hayat/Downloads/R_files/data/significant_enriched_pfam_domains_agelas_combined.tsv", sep='\t', index=False)
# kegg_sig.to_csv(r"C:/Users/hayat/Downloads/R_files/data/enriched_KEGG_Pathways_Agelas_combined.tsv", sep='\t', index=False)

print('Done.')
