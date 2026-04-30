# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 12:11:31 2025

@author: hayat
"""

###############################################
# Scatter Plot: Plasmid Abundance vs Number
# of Plasmids per Sponge Genus, sized by
# Median Plasmid Length
###############################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from matplotlib import rcParams

#import numpy as np
from matplotlib.patches import Patch

# --------------------
# Plot aesthetics
# --------------------
rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.size": 3,
    "xtick.major.width": 0.6,
    "ytick.major.size": 3,
    "ytick.major.width": 0.6,
    "svg.fonttype": "none"
})

# --------------------
# File paths
# --------------------
rpkm_file = "C:/Users/hayat/Downloads/R_files/data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv"
plasmid_names_file = "C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids_names.txt"
host_info_file = "C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt"
plasmid_fasta_file = "C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.fasta"
output_plot_file = "C:/Users/hayat/Downloads/R_files/graphs/rpkm_vs_plasmid_number_per_genus_sized.svg"

# --------------------
# Load data
# --------------------
df = pd.read_csv(rpkm_file, sep="\t", index_col=0)

with open(plasmid_names_file, "r") as f:
    plasmid_names_list = [line.strip() for line in f if line.strip()]

df = df.loc[df.index.isin(plasmid_names_list)]

host_info = pd.read_csv(host_info_file, sep="\t")

# --------------------
# Presence filtering (RPKM ≥ 1)
# --------------------
presence = df >= 1
filtered_df = df[presence.sum(axis=1) >= 0]
filtered_df = filtered_df.loc[:, (filtered_df > 0).any(axis=0)]

# --------------------
# Total RPKM per metagenome
# --------------------
total_rpkm_per_meta = filtered_df.sum(axis=0).reset_index()
total_rpkm_per_meta.columns = ["Metagenome", "Total_RPKM"]

total_rpkm_with_genus = total_rpkm_per_meta.merge(
    host_info, left_on="Metagenome", right_on="Run", how="left"
)
total_rpkm_with_genus = total_rpkm_with_genus[total_rpkm_with_genus["biome_genus"] != "Not_Defined"]

genus_total_rpkm_df = total_rpkm_with_genus.groupby("biome_genus")["Total_RPKM"].sum().reset_index()
genus_total_rpkm_df.columns = ["Sponge Genus", "Total_RPKM"]

# --------------------
# Count plasmids per genus
# --------------------
plasmid_presence = (filtered_df >= 1).astype(int).T

plasmid_with_genus = plasmid_presence.merge(
    host_info, left_index=True, right_on="Run", how="left"
)
plasmid_with_genus = plasmid_with_genus[plasmid_with_genus["biome_genus"] != "Not_Defined"]

plasmids_per_genus = (
    plasmid_with_genus
    .groupby("biome_genus")
    .apply(lambda df: (df.drop(columns=["Run", "biome_genus"]).sum(axis=0) > 0).sum())
    .reset_index()
)
plasmids_per_genus.columns = ["Sponge Genus", "Number of Plasmids"]

# --------------------
# HMA/LMA annotation
# --------------------
hma_lma_dict = {
    'Agelas': 'HMA', 'Aiolochroia': 'HMA', 'Aplysina': 'HMA', 'Chondrilla': 'HMA', 'Coscinoderma': 'HMA',
    'Geodia': 'HMA', 'Ircinia': 'HMA', 'Petrosia': 'HMA', 'Pseudoceratina': 'HMA', 'Rhabdastrella': 'HMA',
    'Sarcotragus': 'HMA', 'Smenospongia': 'HMA', 'Theonella': 'HMA', 'Thoosa': 'HMA', 'Verongula': 'HMA',
    'Rhopaloeides': 'HMA', 'Xestospongia': 'HMA', 'Manihinea': 'HMA',
    'Amphimedon': 'LMA', 'Axinella': 'LMA', 'Baikalospongia': 'LMA', 'Cinachyrella': 'LMA', 'Clathria': 'LMA',
    'Cliona': 'LMA', 'Crella': 'LMA', 'Cymbastela': 'LMA', 'Dysidea': 'LMA', 'Ephydatia': 'LMA', 'Eunapius': 'LMA',
    'Halichondria': 'LMA', 'Haliclona': 'LMA', 'Hymedesmia': 'LMA', 'Ianthella': 'LMA', 'Isodictya': 'LMA',
    'Lamellodysidea': 'LMA', 'Leucetta': 'LMA', 'Mycale': 'LMA', 'Myxilla': 'LMA', 'Niphates': 'LMA',
    'Phyllospongia': 'LMA', 'Scopalina': 'LMA', 'Spheciospongia': 'LMA', 'Spongilla': 'LMA', 'Stylissa': 'LMA',
    'Tedaniidae': 'LMA', 'Pericharax': 'LMA', 'Lophophysema': 'LMA', 'Haplosclerida': 'LMA',
    'Acarnus': 'N.D.', 'Not_Defined': 'N.D.'
}

plasmids_per_genus["HMA_LMA"] = plasmids_per_genus["Sponge Genus"].map(hma_lma_dict)
genus_total_rpkm_df["HMA_LMA"] = genus_total_rpkm_df["Sponge Genus"].map(hma_lma_dict)

# --------------------
# Compute median plasmid length per genus
# --------------------
plasmid_lengths = []
for record in SeqIO.parse(plasmid_fasta_file, "fasta"):
    plasmid_lengths.append({"Plasmid": record.id, "Length": len(record.seq)})

plasmid_lengths_df = pd.DataFrame(plasmid_lengths)

# Melt presence table to long format
plasmid_long = plasmid_with_genus.melt(
    id_vars=["Run", "biome_genus"],
    var_name="Plasmid",
    value_name="Presence"
)
plasmid_long = plasmid_long[plasmid_long["Presence"] >= 1]

# Merge with lengths and compute median
plasmid_long = plasmid_long.merge(plasmid_lengths_df, on="Plasmid", how="left")
median_length_per_genus = plasmid_long.groupby("biome_genus")["Length"].median().reset_index()
median_length_per_genus.columns = ["Sponge Genus", "Median_Plasmid_Length"]

# --------------------
# Prepare final scatter plot DataFrame
# --------------------
scatter_df = plasmids_per_genus.merge(
    genus_total_rpkm_df[["Sponge Genus", "Total_RPKM"]],
    on="Sponge Genus",
    how="inner"
)
scatter_df = scatter_df.merge(median_length_per_genus, on="Sponge Genus", how="left")

# Only HMA/LMA
scatter_df = scatter_df[scatter_df["HMA_LMA"].isin(["HMA", "LMA"])]

# Split by HMA/LMA
hma_df = scatter_df[scatter_df["HMA_LMA"] == "HMA"]
lma_df = scatter_df[scatter_df["HMA_LMA"] == "LMA"]

# Select top 5 in each group (by Number of Plasmids)
top_hma = hma_df.nlargest(5, "Number of Plasmids")
top_lma = lma_df.nlargest(5, "Number of Plasmids")
top_labels = pd.concat([top_hma, top_lma])


plt.figure(figsize=(8,6))

# Plot HMA first
sns.scatterplot(
    data=scatter_df[scatter_df["HMA_LMA"]=="HMA"],
    x="Total_RPKM",
    y="Number of Plasmids",
    size="Median_Plasmid_Length",
    sizes=(50, 500),
    color="darkgreen",
    edgecolor="k",
    alpha=0.8,
    zorder=1
)

# Plot LMA second (on top)
sns.scatterplot(
    data=scatter_df[scatter_df["HMA_LMA"]=="LMA"],
    x="Total_RPKM",
    y="Number of Plasmids",
    size="Median_Plasmid_Length",
    sizes=(50, 500),
    color="orange",
    edgecolor="k",
    alpha=0.8,
    zorder=2
)

# Add labels for top 5 HMA + top 5 LMA
for i, row in top_labels.iterrows():
    plt.text(
        row["Total_RPKM"], row["Number of Plasmids"], row["Sponge Genus"],
        fontsize=9, weight="bold"
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Total Plasmid RPKM per Genus (log scale)")
plt.ylabel("Number of Plasmids per Genus (log scale)")
plt.grid(True, which="both", linestyle="--", alpha=0.3)

# ---- Manual legend ----
# HMA/LMA color legend
color_legend = [
    Patch(facecolor="darkgreen", edgecolor="k", label="HMA"),
    Patch(facecolor="orange", edgecolor="k", label="LMA")
]

# Compute circle sizes exactly as in the plot
min_size, max_size = 50, 500
lengths = scatter_df["Median_Plasmid_Length"]
scatter_df["Plot_Size"] = min_size + (lengths - lengths.min()) / (lengths.max() - lengths.min()) * (max_size - min_size)

# Pick representative lengths for size legend (min, median, max)
size_legend_values = [lengths.min(), lengths.median(), lengths.max()]
size_legend_handles = [
    plt.scatter([], [], s=min_size + (v - lengths.min()) / (lengths.max() - lengths.min()) * (max_size - min_size),
                color="gray", edgecolor="k") for v in size_legend_values
]

# Add legends
first_legend = plt.legend(handles=color_legend, title="Host Type", loc="upper left")
plt.gca().add_artist(first_legend)
plt.legend(handles=size_legend_handles,
           labels=[f"{int(v)} bp" for v in size_legend_values],
           title="Median Plasmid Length",
           loc="upper right",
           scatterpoints=1)


plt.tight_layout()
plt.savefig(output_plot_file)
plt.show()






# Save scatter plot data to CSV
scatter_df.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/scatter_plot_data_length_RPKM_plasmid_number_per_genus.csv",
    index=False
)
