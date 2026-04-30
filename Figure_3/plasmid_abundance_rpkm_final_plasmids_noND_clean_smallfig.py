###############################################
# Plasmid RPKM Analysis Pipeline
#
# Purpose:
#   - Filter plasmid RPKM data
#   - Summarize plasmid abundance and presence
#   - Generate plots of abundance across sponge genera
#   - Identify widespread & abundant plasmids
#   - Compare plasmid abundance/diversity between HMA vs LMA sponges
#
# Inputs:
#   - CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv
#   - all_sponge_plasmids_names.txt
#   - Accessions_to_SpongeGenusIDs.txt
#
# Outputs (commented out, can be enabled as needed):
#   - CSV exports of filtered and summary tables
#   - Plots (PNG/SVG)
###############################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#import numpy as np
from matplotlib.patches import Patch
from collections import defaultdict
from scipy.stats import mannwhitneyu
from matplotlib import rcParams
#from svgutils.compose import Figure, Panel, SVG, Text
#from adjustText import adjust_text
rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.size": 3,
    "xtick.major.width": 0.6,
    "ytick.major.size": 3,
    "ytick.major.width": 0.6,
    "svg.fonttype": "none"   # <-- keeps text as editable text
})

# --------------------
# Load & filter data
# --------------------

df = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv",
    sep="\t",
    index_col=0
)

with open("C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids_names.txt", "r") as f:
    plasmid_names_list = [line.strip() for line in f if line.strip()]

# keep only plasmids of interest
df = df.loc[df.index.isin(plasmid_names_list)]

# --------------------
# Presence filtering
# --------------------

presence = df >= 1
filtered_df = df[presence.sum(axis=1) >= 0]
filtered_df = filtered_df.loc[:, (filtered_df > 0).any(axis=0)]  # drop all-zero metagenomes

# --------------------
# Abundance summaries
# --------------------

most_abundant_by_sample = filtered_df.idxmax()
abundance_values = filtered_df.max()

abundance_summary = pd.DataFrame({
    "most_abundant_plasmid": most_abundant_by_sample,
    "RPKM": abundance_values
}).reset_index().rename(columns={"index": "Metagenome"})

# --------------------
# Merge with host info
# --------------------

host_info = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt",
    sep="\t"
)



merged_df = abundance_summary.merge(host_info, left_on="Metagenome", right_on="Run", how="left")
merged_df = merged_df.drop(columns=["Run"])


# merged_df.to_csv(".../abundance_with_host.csv", index=False)

print("Merged abundance summary:")
print(merged_df.head())

# --------------------
# HMA/LMA dictionary and colors
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

# --------------------
# Genus-level summaries
# --------------------

# Count metagenomes per genus
genus_counts = host_info["biome_genus"].value_counts().reset_index()
genus_counts.columns = ["Sponge Genus", "Number of Metagenomes"]
genus_counts = genus_counts.sort_values("Number of Metagenomes", ascending=False)

plt.figure(figsize=(10, 6))
sns.barplot(data=genus_counts, x="Sponge Genus", y="Number of Metagenomes", palette="viridis")
plt.title("Number of Metagenomes per Sponge Genus")
plt.ylabel("Number of Metagenomes")
plt.xlabel("Sponge Genus")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
# plt.savefig(".../metagenome_count_per_genus.png", dpi=300)
plt.show()

# Count plasmids per genus

plasmid_presence = (filtered_df >= 1).astype(int).T

# Merge but keep all rows for now
plasmid_with_genus = plasmid_presence.merge(
    host_info, left_index=True, right_on="Run", how="left"
)

# ---- NEW LINE: filter out Not_Defined ONLY FOR THIS PLOT ----
plasmid_with_genus = plasmid_with_genus[plasmid_with_genus["biome_genus"] != "Not_Defined"]

# Compute plasmids per genus
plasmids_per_genus = (
    plasmid_with_genus
    .groupby("biome_genus")
    .apply(lambda df: (df.drop(columns=["Run", "biome_genus"]).sum(axis=0) > 0).sum())
    .reset_index()
)

plasmids_per_genus.columns = ["Sponge Genus", "Number of Plasmids"]

# Add HMA/LMA annotation
plasmids_per_genus["HMA_LMA"] = plasmids_per_genus["Sponge Genus"].map(hma_lma_dict)

# Sort
plasmids_per_genus = plasmids_per_genus.sort_values("Number of Plasmids", ascending=False)

# Plot
plt.figure(figsize=(10, 6))
sns.barplot(
    data=plasmids_per_genus,
    x="Sponge Genus", y="Number of Plasmids",
    hue="HMA_LMA",
    palette={"HMA": "darkgreen", "LMA": "orange", "N.D.": "lightgrey"}
)
plt.title("Number of Plasmids per Sponge Genus (RPKM ≥ 1)")
plt.ylabel("Number of Plasmids")
plt.xlabel("Sponge Genus")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_count_per_genus.svg")
plt.show()

# Save the plasmids per genus table as CSV
plasmids_per_genus.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_count_per_genus.csv",
    index=False
)



# # --------------------
# # Abundance per genus (wrong!!!! Redo)
# # --------------------

# genus_total_rpkm = merged_df.groupby("biome_genus")["RPKM"].sum().sort_values(ascending=False)

# log_total_rpkm = np.log10(genus_total_rpkm)

# plt.figure(figsize=(12, 6))
# sns.barplot(x=log_total_rpkm.index, y=genus_total_rpkm.values, palette="muted")
# plt.xticks(rotation=45, ha="right")
# plt.ylabel("Total RPKM")
# plt.title("Total Plasmid Abundance per Sponge Genus")
# plt.xlabel("Sponge Genus")
# plt.tight_layout()
# # plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/rpkm_count_per_genus.png", dpi=300)
# plt.show()

# --------------------
# TRUE total RPKM per metagenome
# --------------------

# Sum all plasmid RPKMs per metagenome
total_rpkm_per_meta = filtered_df.sum(axis=0).reset_index()
total_rpkm_per_meta.columns = ["Metagenome", "Total_RPKM"]

# Merge with host_info to get sponge genus/species
total_rpkm_with_genus = total_rpkm_per_meta.merge(
    host_info, left_on="Metagenome", right_on="Run", how="left"
)

# ---- NEW LINE: remove Not_Defined ONLY for this genus plot ----
total_rpkm_with_genus = total_rpkm_with_genus[
    total_rpkm_with_genus["biome_genus"] != "Not_Defined"
]

# --------------------
# TRUE total RPKM per sponge genus
# --------------------
genus_total_rpkm_true = (
    total_rpkm_with_genus.groupby("biome_genus")["Total_RPKM"]
    .sum()
    .sort_values(ascending=False)
)

genus_total_rpkm_df = genus_total_rpkm_true.reset_index()
genus_total_rpkm_df.columns = ["Sponge Genus", "Total_RPKM"]
genus_total_rpkm_df["HMA_LMA"] = genus_total_rpkm_df["Sponge Genus"].map(hma_lma_dict)

plt.figure(figsize=(12, 6))
sns.barplot(
    data=genus_total_rpkm_df,
    x="Sponge Genus", y="Total_RPKM",
    hue="HMA_LMA",
    palette={"HMA": "darkgreen", "LMA": "orange", "N.D.": "lightgrey"}
)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Total RPKM")
plt.title("Total Plasmid Abundance per Sponge Genus")
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/rpkm_count_per_genus.svg")
plt.show()

# Save total RPKM per genus as CSV
genus_total_rpkm_df.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/rpkm_count_per_genus.csv",
    index=False
)

# --------------------
# Prepare data for scatter plot
# --------------------
# Merge plasmid count per genus with total RPKM per genus
scatter_df = plasmids_per_genus.merge(
    genus_total_rpkm_df[["Sponge Genus", "Total_RPKM"]],
    on="Sponge Genus",
    how="inner"
)


# Exclude Not_Defined just to be safe
scatter_df = scatter_df[scatter_df["HMA_LMA"].isin(["HMA", "LMA"])]

# --------------------
# Scatter plot
# --------------------
plt.figure(figsize=(8,6))
sns.scatterplot(
    data=scatter_df,
    x="Total_RPKM",
    y="Number of Plasmids",
    hue="HMA_LMA",
    palette={"HMA": "darkgreen", "LMA": "orange"},
    s=100,
    edgecolor="k",
    alpha=0.8
)

texts = []
for i, row in scatter_df.iterrows():
    texts.append(
        plt.text(row["Total_RPKM"], row["Number of Plasmids"], row["Sponge Genus"], fontsize=9)
    )

#adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Total Plasmid RPKM per Genus (log scale)")
plt.ylabel("Total Number of Plasmids per Genus (log scale)")
#plt.title("Plasmid Number vs Abundance per Host Genus")
plt.legend(title="Host Type")
plt.grid(True, which="both", linestyle="--", alpha=0.3)
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/rpkm_vs_plasmid_number_per_genus.svg")
plt.show()




# --------------------
# Load mapping
# --------------------
mapping_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_id_mapping.csv")
mapping_dict = mapping_df.set_index("Original_ID")["New_ID"].to_dict()

# --------------------
# Widespread & abundant plasmids
# --------------------

min_presence = int(0.01 * df.shape[1])
filtered_df_1 = df[presence.sum(axis=1) >= min_presence]

presence_counts = (filtered_df_1 >= 1).sum(axis=1)
total_abundance = filtered_df_1.sum(axis=1)

presence_thresh = presence_counts.quantile(0.90)
total_abundance_thresh = total_abundance.quantile(0.90)

most_widespread = presence_counts[presence_counts >= presence_thresh]
top_abundant_total = total_abundance[total_abundance >= total_abundance_thresh]

top_both = set(most_widespread.index) & set(top_abundant_total.index)
print(f"{len(top_both)} plasmids are both highly widespread and abundant (top 10%).")


# --------------------
# Identify top 5 widespread + top 5 abundant
# --------------------

top5_widespread = presence_counts.sort_values(ascending=False).head(5).index
top5_abundant = total_abundance.sort_values(ascending=False).head(5).index

labels_to_annotate = set(top5_widespread) | set(top5_abundant)


# --------------------
# Scatterplot
# --------------------
plt.figure(figsize=(3.6, 2.5), dpi=96)
plt.scatter(
    presence_counts,
    total_abundance + 1e-3,
    c=["blue" if idx in top_both else "gray" for idx in presence_counts.index],
    s=[10 if idx in top_both else 3 for idx in presence_counts.index],
    alpha=0.7, edgecolor='k', linewidth=0.3
)

plt.axhline(y=total_abundance_thresh, color='darkgray', linestyle='--')
plt.axvline(x=presence_thresh, color='darkblue', linestyle='--')
plt.yscale("log")
plt.xlabel("Number of Metagenomes with RPKM ≥ 1")
plt.ylabel("Total RPKM across Metagenomes (log scale)")
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.3)
plt.tight_layout()

# --------------------
# Annotate top 5 widespread & top 5 abundant
# --------------------
for plasmid in labels_to_annotate:
    x = presence_counts.loc[plasmid]
    y = total_abundance.loc[plasmid]

    label = mapping_dict.get(plasmid, plasmid)  # Use New_ID if exists

    plt.text(
        x + 0.3, y * 1.05,       # small offset
        label,
        fontsize=6,
        color="black"
    )
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_vs_abundance_highlighted_final_plasmids_annotated.svg")
plt.show()


# --------------------
# Summary table
# --------------------
widespread_abundant_df = filtered_df_1.loc[list(top_both)]

summary_df = pd.DataFrame({
    "Original_ID": widespread_abundant_df.index,
    "Presence_Count": (widespread_abundant_df >= 1).sum(axis=1),
    "Mean_RPKM": widespread_abundant_df.mean(axis=1),
    "Total_RPKM": widespread_abundant_df.sum(axis=1)
}).reset_index(drop=True)

summary_df = summary_df.merge(mapping_df, on="Original_ID", how="left")

summary_df = summary_df[["Original_ID", "New_ID", "Presence_Count", "Total_RPKM", "Mean_RPKM"]]

#summary_df.to_csv(
#    "C:/Users/hayat/Downloads/R_files/graphs/widespread_abundant_plasmids_summary.csv",
#    index=False
#)





# --------------------
# Genus filtering for clustermap
# --------------------

top_both_df = filtered_df_1.loc[filtered_df.index.intersection(top_both)]

# map genus
genus_map = merged_df[["Metagenome", "biome_genus"]].copy()
genus_map = genus_map.rename(columns={"biome_genus": "Genus"})
genus_counts = genus_map["Genus"].value_counts()

valid_genera = genus_counts[genus_counts > 5].index
genus_map = genus_map[genus_map["Genus"].isin(valid_genera)]
genus_map["genus_with_count"] = genus_map["Genus"] + " (n=" + genus_map["Genus"].map(genus_counts).astype(str) + ")"

merged_df = merged_df.drop(columns=["genus_with_count"], errors="ignore")
merged_df = merged_df.merge(genus_map[["Metagenome", "genus_with_count"]], on="Metagenome", how="left")

# prepare long format for clustermap
long_df = top_both_df.reset_index().melt(id_vars="Genome", var_name="Metagenome", value_name="RPKM")
long_df = long_df.rename(columns={"Genome": "Plasmid"})
long_df = long_df.merge(merged_df[["Metagenome", "genus_with_count"]], on="Metagenome", how="left")
long_df = long_df[long_df["RPKM"] > 0]

# unique plasmid-metagenome combos
long_df_unique = long_df.drop_duplicates(subset=['Plasmid', 'Metagenome'])
presence_count_matrix = long_df_unique.pivot_table(index='Plasmid', columns='genus_with_count', aggfunc='size', fill_value=0)

metagenomes_per_genus = merged_df[["Metagenome", "genus_with_count"]].drop_duplicates()["genus_with_count"].value_counts()
metagenomes_per_genus = metagenomes_per_genus.reindex(presence_count_matrix.columns)

relative_presence_matrix = presence_count_matrix.divide(metagenomes_per_genus, axis=1) * 100



hma_lma_colors = {"HMA": "darkgreen", "LMA": "orange", "N.D.": "lightgrey"}

filtered_columns = [col for col in relative_presence_matrix.columns if hma_lma_dict.get(col.split(" (n=")[0]) != "N.D."]
relative_presence_matrix = relative_presence_matrix[filtered_columns]

col_colors = [hma_lma_colors.get(hma_lma_dict.get(col.split(" (n=")[0], "N.D."), "lightgrey") for col in filtered_columns]

# --------------------
# Clustermap plot
# --------------------

# relative_presence_matrix.to_csv(".../top_abundant_widespread_plasmid_relative_presence_matrix_filtered.csv")

g = sns.clustermap(
    relative_presence_matrix,
    cmap="Blues",
    figsize=(6.5, 4.3),
    cbar_kws={"label": "Presence (% of metagenomes)"},
    method="average",
    metric="euclidean",
    dendrogram_ratio=(0.1, 0.1),
    colors_ratio=0.01,
    col_colors=col_colors
)

#plt.suptitle("Clustered Relative Presence of Top 10% Abundant & Widespread Plasmids\nAcross Host Genera (>1 metagenome with plasmid)", fontsize=18, y=1.05)

legend_handles = [Patch(color=color, label=label) for label, color in hma_lma_colors.items() if label != "N.D."]
g.ax_col_dendrogram.legend(handles=legend_handles, title="Host Microbial Abundance", loc="center", bbox_to_anchor=(1.2, -3), ncol=3)

g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="center")
g.ax_heatmap.tick_params(axis='y', which='both', length=0)
g.ax_heatmap.set_yticklabels([])
g.ax_heatmap.set_xlabel("Host Genus (with metagenome count)")
g.ax_heatmap.set_ylabel("Plasmid")

g.cax.set_position([0.93, 0.2, 0.02, 0.3])  # [x, y, width, height]
g.cax.yaxis.set_label_position("right")
g.cax.yaxis.tick_right()

# g.savefig(".../top_10_percent_abundant_widespread_relative_presence_clustermap.png", dpi=300)
#g.savefig("C:/Users/hayat/Downloads/R_files/graphs/publication_figures/top_10_percent_abundant_widespread_relative_presence_clustermap_Figure_2d.svg")
plt.show()

# --------------------
# Plasmids vs sponge genera matrix
# --------------------
# Create plasmid presence/absence matrix (1 if RPKM > 0, else 0)
#plasmid_presence = (filtered_df > 0).astype(int)

metagenome_to_genus = host_info.set_index("Run")["biome_genus"].to_dict()
genus_to_metagenomes = defaultdict(list)
for metagenome, genus in metagenome_to_genus.items():
    if genus:
        genus_to_metagenomes[genus].append(metagenome)

plasmid_vs_genus = pd.DataFrame(0, index=plasmid_presence.index, columns=genus_to_metagenomes.keys())
for genus, metagenomes in genus_to_metagenomes.items():
    valid_metagenomes = [m for m in metagenomes if m in plasmid_presence.columns]
    genus_presence = (plasmid_presence[valid_metagenomes].sum(axis=1) > 0).astype(int)
    plasmid_vs_genus[genus] = genus_presence

#plasmid_vs_genus.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmids_vs_sponge_genera.csv")

# --------------------
# RPKM per HMA vs LMA
# --------------------

rpkm_per_metagenome = filtered_df.sum(axis=0).reset_index()
rpkm_per_metagenome.columns = ["Metagenome", "Total_RPKM"]

rpkm_with_genus = rpkm_per_metagenome.merge(host_info, left_on="Metagenome", right_on="Run", how="left")
rpkm_with_genus["Type"] = rpkm_with_genus["biome_genus"].map(hma_lma_dict)
rpkm_with_genus = rpkm_with_genus.dropna(subset=["Type"])
rpkm_with_genus = rpkm_with_genus[rpkm_with_genus["Type"].isin(["HMA", "LMA"])]

order = ["HMA", "LMA"]
hma_data = rpkm_with_genus[rpkm_with_genus["Type"] == "HMA"]["Total_RPKM"]
lma_data = rpkm_with_genus[rpkm_with_genus["Type"] == "LMA"]["Total_RPKM"]
stat, p_value = mannwhitneyu(hma_data, lma_data, alternative="two-sided")

plt.figure(figsize=(2.2, 1.8), dpi=96)
sns.boxplot(data=rpkm_with_genus, x="Type", y="Total_RPKM", palette=hma_lma_colors, showfliers=True, order=order)

means = rpkm_with_genus.groupby("Type")["Total_RPKM"].mean()
medians = rpkm_with_genus.groupby("Type")["Total_RPKM"].median()

summary_rpkm = pd.DataFrame({
    "Type": order,
    "Mean_Total_RPKM": [means[g] for g in order],
    "Median_Total_RPKM": [medians[g] for g in order],
    "P_Value_Total_RPKM": [p_value] + [None]
})

for i, g in enumerate(order):
    plt.scatter(i, means[g], color='blue', marker='D', s=15, edgecolor='black', label='Mean' if i == 0 else "", zorder=10)
    plt.scatter(i, medians[g], color='red', marker='o', s=15, edgecolor='black', label='Median' if i == 0 else "", zorder=10)

plt.yscale("log")
#plt.title("Total Plasmid RPKM per Metagenome by Sponge Type (HMA vs LMA)")
plt.ylabel("Total RPKM (log scale)")
plt.xlabel("Sponge Type")
plt.legend()

# p-value annotation
x1, x2 = 0, 1
y = rpkm_with_genus["Total_RPKM"].max() * 1.1
plt.plot([x1, x1, x2, x2], [y, y+0.2, y+0.2, y], lw=1.5, c='k')
plt.text(0.5, y+0.25, "****" if p_value < 0.0001 else "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns", ha='center')

# plt.savefig(".../plasmid_RPKM_per_hma_lma.png", dpi=300)
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_RPKM_per_hma_lma_Figure_3b.svg")
plt.show()

# --------------------
# Plasmid richness per HMA vs LMA
# --------------------

plasmid_richness = (filtered_df >= 1).sum(axis=0).reset_index()
plasmid_richness.columns = ["Metagenome", "Plasmid_Richness"]

richness_with_genus = plasmid_richness.merge(host_info, left_on="Metagenome", right_on="Run", how="left")
richness_with_genus["Type"] = richness_with_genus["biome_genus"].map(hma_lma_dict)
richness_with_genus = richness_with_genus.dropna(subset=["Type"])
richness_with_genus = richness_with_genus[richness_with_genus["Type"].isin(["HMA", "LMA"])]

hma_rich = richness_with_genus[richness_with_genus["Type"] == "HMA"]["Plasmid_Richness"]
lma_rich = richness_with_genus[richness_with_genus["Type"] == "LMA"]["Plasmid_Richness"]
stat_rich, p_rich = mannwhitneyu(hma_rich, lma_rich, alternative="two-sided")

plt.figure(figsize=(2.2, 1.8), dpi=96)

# Boxplot
sns.boxplot(
    data=richness_with_genus,
    x="Type",
    y="Plasmid_Richness",
    palette=hma_lma_colors,
    showfliers=True,
    order=order
)

# Calculate means and medians
grouped = richness_with_genus.groupby("Type")["Plasmid_Richness"]
means = grouped.mean()
medians = grouped.median()

# Overlay mean and median markers
for i, g in enumerate(order):
    plt.scatter(
        i, means[g],
        color='blue',
        marker='D',
        label='Mean' if i == 0 else "",
        s=15,
        edgecolor='black',
        zorder=10
    )
    plt.scatter(
        i, medians[g],
        color='red',
        marker='o',
        label='Median' if i == 0 else "",
        s=15,
        edgecolor='black',
        zorder=10
    )

#plt.title("Plasmid Diversity per Metagenome by Sponge Type (HMA vs LMA)")
plt.ylabel("Plasmid Richness (RPKM ≥ 1)")
plt.xlabel("Sponge Type")
plt.legend()

# Annotate p-value
x1, x2 = 0, 1
y = max(hma_rich.max(), lma_rich.max()) * 1.05
h = 2
plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c='k')

if p_rich < 0.0001:
    p_text = "****"
elif p_rich < 0.001:
    p_text = "***"
elif p_rich < 0.01:
    p_text = "**"
elif p_rich < 0.05:
    p_text = "*"
else:
    p_text = "ns"

plt.text((x1 + x2) * 0.5, y + h + 1, p_text, ha='center', va='bottom', color='k')

#plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/publication_figures/plasmid_diversity_per_hma_lma_figure_3a.svg")
plt.show()


means = richness_with_genus.groupby("Type")["Plasmid_Richness"].mean()
medians = richness_with_genus.groupby("Type")["Plasmid_Richness"].median()

summary_richness = pd.DataFrame({
    "Type": order,
    "Mean_Plasmid_Richness": [means[g] for g in order],
    "Median_Plasmid_Richness": [medians[g] for g in order],
    "P_Value_Plasmid_Richness": [p_rich] + [None]*(len(order)-1)

})




# Merge both summaries on "Type"
summary_all = pd.merge(summary_rpkm, summary_richness, on="Type")

# Save as CSV
#summary_all.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_rpkm_richness_summary.csv", index=False)

# Optional: save as TSV
#summary_all.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_rpkm_richness_summary.tsv", sep="\t", index=False)

# ---- Outliers_abundance ----

lma_df = rpkm_with_genus[rpkm_with_genus["Type"] == "LMA"]

# ---- 2. Compute IQR for Total_RPKM ----
Q1 = lma_df["Total_RPKM"].quantile(0.25)
Q3 = lma_df["Total_RPKM"].quantile(0.75)
IQR = Q3 - Q1

lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

# ---- 3. Select outliers ----
outliers = lma_df[(lma_df["Total_RPKM"] < lower_bound) | 
                  (lma_df["Total_RPKM"] > upper_bound)]

# ---- 4. Print results ----
print("Number of outlier metagenomes in LMA:", outliers.shape[0])
print("\nOutlier metagenomes with species:")
print(outliers[["Metagenome", "Total_RPKM", "biome_genus"]])


# --------------------
# Detect outliers in LMA plasmid richness
# --------------------

# Filter to LMA only
lma_df = richness_with_genus[richness_with_genus["Type"] == "LMA"].copy()

# Compute IQR outlier thresholds
Q1 = lma_df["Plasmid_Richness"].quantile(0.25)
Q3 = lma_df["Plasmid_Richness"].quantile(0.75)
IQR = Q3 - Q1

lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

# Identify low or high outliers
lma_outliers = lma_df[
    (lma_df["Plasmid_Richness"] < lower_bound) |
    (lma_df["Plasmid_Richness"] > upper_bound)
].copy()

print("Number of LMA outliers:", lma_outliers.shape[0])

# Show which metagenomes and species they belong to
print("\nLMA Outlier Metagenomes + Host Species:")
print(lma_outliers)
