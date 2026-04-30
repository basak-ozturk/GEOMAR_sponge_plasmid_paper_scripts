# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 09:11:25 2025

@author: hayat
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.colors as mcolors
#from matplotlib.cm import get_cmap
import seaborn as sns

plt.rcParams["svg.fonttype"] = "none"  # keep text as text in SVG
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 8
sns.set_style("whitegrid")  # seaborn white background
sns.set_context("notebook", font_scale=1.0)  # keep 8 pt


# === 1. Load the data ===
cyto_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/cytoscape_edges_with_mcl_clusters.csv")
nodes_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/cytoscape_nodes.csv")


# === 2. Normalize cluster names, extract singletons ===
cyto_df['Cluster'] = cyto_df['__mclCluster'].astype(str)
nodes_df['name'] = nodes_df['ID']
cyto_plasmids = set(cyto_df['name'])
nodes_df['is_singleton'] = ~nodes_df['name'].isin(cyto_plasmids)
singleton_df = nodes_df[nodes_df['is_singleton']].copy()
singleton_df['Cluster'] = 'Singleton'

# === 3. Harmonize and combine ===
# Normalize case in Mob and Integrase columns (lowercase: 'yes'/'no')
cyto_df['Mob'] = cyto_df['Mob'].str.lower()
cyto_df['Integrase'] = cyto_df['Integrase'].str.lower()
nodes_df['Mob'] = nodes_df['Mob'].str.lower()
nodes_df['Integrase'] = nodes_df['Integrase'].str.lower()

singleton_df = singleton_df[['Cluster', 'Mob', 'Integrase', 'name']]
cyto_df_sub = cyto_df[['Cluster', 'Mob', 'Integrase', 'name']]
combined_df = pd.concat([cyto_df_sub, singleton_df], ignore_index=True)

# === 4. Filter clusters with ≥5 plasmids ===
cluster_sizes = combined_df['Cluster'].value_counts()
valid_clusters = cluster_sizes[cluster_sizes >= 5].index
filtered_df = combined_df[combined_df['Cluster'].isin(valid_clusters)].copy()

# === 5. Fix missing categories using reindexing ===
# Mob counts
mob_counts = (
    filtered_df.groupby(['Cluster', 'Mob'])
    .size()
    .unstack()
    .rename(columns={'yes': 'Mob+', 'no': 'Mob-'})
    .reindex(columns=['Mob+', 'Mob-'], fill_value=0)
)

# Integrase counts
int_counts = (
    filtered_df.groupby(['Cluster', 'Integrase'])
    .size()
    .unstack()
    .rename(columns={'yes': 'Int+', 'no': 'Int-'})
    .reindex(columns=['Int+', 'Int-'], fill_value=0)
)

# === 6. Order clusters: descending size, singleton last ===
cluster_order = mob_counts.sum(axis=1).sort_values(ascending=False).index.tolist()
if "Singleton" in cluster_order:
    cluster_order = [c for c in cluster_order if c != "Singleton"] + ["Singleton"]
mob_counts = mob_counts.loc[cluster_order]
int_counts = int_counts.loc[cluster_order]

# === 7. Prepare for horizontal barplot with Seaborn ===
# We "melt" both tables into long format, and label them as Mob/Integrase

mob_long = mob_counts.reset_index().melt(id_vars="Cluster", var_name="Category", value_name="Count")
mob_long["Type"] = "Mob"

int_long = int_counts.reset_index().melt(id_vars="Cluster", var_name="Category", value_name="Count")
int_long["Type"] = "Integrase"

plot_df = pd.concat([mob_long, int_long], ignore_index=True)

# Combine Type and Category for better coloring
plot_df["CatLabel"] = plot_df["Type"] + ": " + plot_df["Category"]

# Offset bars to group by cluster with spacing between Mob and Integrase
plot_df["Cluster_Type"] = plot_df["Cluster"] + " (" + plot_df["Type"] + ")"

# === 8. Plot ===
#ns.set(style="whitegrid", font_scale=1.2)
# --- Color palette and ordering ---
palette = sns.color_palette("viridis", n_colors=4)
cat_order = ["Mob: Mob-", "Mob: Mob+", "Integrase: Int-", "Integrase: Int+"]

plot_df["CatLabel"] = pd.Categorical(plot_df["CatLabel"], categories=cat_order, ordered=True)

# Prepare y order: Mob first, Integrase second per cluster
cluster_order_expanded = []
for cl in cluster_order:
    cluster_order_expanded.append(f"{cl} (Mob)")
    cluster_order_expanded.append(f"{cl} (Integrase)")

plot_df["Cluster_Type"] = pd.Categorical(plot_df["Cluster_Type"], categories=cluster_order_expanded, ordered=True)

# Pivot for stacked bar
pivot_df = plot_df.pivot_table(index="Cluster_Type", columns="CatLabel", values="Count", fill_value=0)
pivot_df.to_csv("C:/Users/hayat/Downloads/R_files/data/mob_integrase_cluster_counts_cytoscape.csv")

# --- Figure size ---
fig_height = max(2, 0.3 * len(pivot_df))  # 0.3 in per row, min 2 in
fig, ax = plt.subplots(figsize=(4, fig_height))  # 100 mm width ~ 4 in

# --- Plot horizontal stacked bars ---
bottom = np.zeros(len(pivot_df))
for i, cat in enumerate(cat_order):
    ax.barh(pivot_df.index, pivot_df[cat], left=bottom, color=palette[i], label=cat)
    bottom += pivot_df[cat].values

# --- Labels ---
ax.set_xlabel("Number of Plasmids", fontsize=8)
ax.set_ylabel("Cluster", fontsize=8)
ax.tick_params(axis='y', labelsize=8)
ax.tick_params(axis='x', labelsize=8)

# --- Legend ---
ax.legend(title="Category", title_fontsize=8, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")

# --- Layout and save ---
plt.tight_layout()
plt.savefig(
    "C:/Users/hayat/Downloads/R_files/graphs/mob_integrase_cluster_distribution.svg",
    format="svg",
    bbox_inches="tight"
)
plt.show()