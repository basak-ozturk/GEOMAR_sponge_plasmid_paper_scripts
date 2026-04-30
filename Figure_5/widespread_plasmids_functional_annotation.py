import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.size": 3,
    "xtick.major.width": 0.6,
    "ytick.major.size": 3,
    "ytick.major.width": 0.6,
})
plt.rcParams['svg.fonttype'] = 'none'

# Load pre-filtered annotation file
annotation_file = "C:/Users/hayat/Downloads/R_files/data/filtered_foreground_eggnog_annotations.tsv"

# Read into DataFrame, skipping comment lines
df = pd.read_csv(annotation_file, sep='\t')

# Fill missing values
df["PFAMs"] = df["PFAMs"].fillna("").astype(str)
df["Description"] = df["Description"].fillna("").astype(str)

# Align PFAMs with Descriptions
pfam_entries = []
for _, row in df.iterrows():
    pfams = [pfam.strip() for pfam in row["PFAMs"].split(",") if pfam.strip() and pfam.strip() != "-"]
    for pfam in pfams:
        pfam_entries.append((pfam, row["Description"]))

pfam_desc_df = pd.DataFrame(pfam_entries, columns=["PFAM", "Description"])
pfam_counts = pfam_desc_df["PFAM"].value_counts().reset_index()
pfam_counts.columns = ["PFAM", "Count"]

# Top PFAMs
top_n = 20
top_pfam_df = pfam_counts.head(top_n)
#top_pfam_df.to_csv("C:/Users/hayat/Downloads/R_files/data/top_pfam_for_curation.csv", index=False)

# Load curated PFAM categories
top_pfam_df_curated = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_pfam_for_curation_curated_consolidated.csv")
top_pfam_df_curated = top_pfam_df_curated.sort_values("Count", ascending=False)
top_20_pfam_df = top_pfam_df_curated.head(20)
category_palette = {
    "Toxin-antitoxin": "#D55E00",               # reddish-orange
    "Oxidoreductase": "#0072B2",                # blue
    "Signalling": "#009E73",                    # bluish green
    "Broad/uncharacterized function": "#999999",# gray
    "Restriction-DNA modification": "#E69F00",  # orange
    "Conjugation": "#56B4E9",                   # light blue
    "Transposase": "#F0E442",                   # yellow
    "Integrase": "#CC79A7",                     # magenta
    "Regulatory protein": "#882255",            # deep burgundy
    "Transporter": "#44AA99",                   # teal
    "Secretion system": "#117733"               # dark green
}



# PFAM plot
plt.figure(figsize=(4, 3))
sns.set_style("whitegrid")

ax = sns.barplot(
    data=top_20_pfam_df,
    y="PFAM", x="Count",
    hue="Category",
    dodge=False,
    palette=category_palette

)

#plt.title("Top PFAM Domains in Widespread Plasmid Proteins (Categorized)", fontsize=16)
plt.xlabel("Count")
plt.ylabel("PFAM")
#plt.yticks(fontsize=8)
#plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.legend(title="Category", loc='upper right', fontsize=8, title_fontsize=8, frameon=False)

plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_plasmids_top_pfam_plot_categorized.png", dpi=300)
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_plasmids_top_pfam_plot_categorized.svg")

plt.show()

# COG processing
cog_series = df["COG_category"].dropna().astype(str).apply(list).explode()
cog_series = cog_series[~cog_series.isin(["-", "S"])]  # Exclude unwanted

cog_counts = cog_series.value_counts().reset_index()
cog_counts.columns = ["COG", "Count"]

# Map COG to functional names
cog_function_map = {
    "J": "Translation (J)",
    "A": "RNA processing and modification (A)",
    "K": "Transcription (K)",
    "L": "Replication, recombination and repair (L)",
    "B": "Chromatin structure and dynamics (B)",
    "D": "Cell cycle control (D)",
    "Y": "Nuclear structure (Y)",
    "V": "Defense mechanisms (V)",
    "T": "Signal transduction (T)",
    "M": "Cell wall/membrane biogenesis (M)",
    "N": "Cell motility (N)",
    "U": "Intracellular trafficking (U)",
    "O": "Posttranslational modification (O)",
    "C": "Energy production (C)",
    "G": "Carbohydrate metabolism (G)",
    "E": "Amino acid metabolism (E)",
    "F": "Nucleotide metabolism (F)",
    "H": "Coenzyme metabolism (H)",
    "I": "Lipid metabolism (I)",
    "P": "Inorganic ion transport (P)",
    "Q": "Secondary metabolites biosynthesis (Q)",
    "R": "General function prediction (R)",
}

cog_counts["Function"] = cog_counts["COG"].map(cog_function_map)

# COG plot
plt.figure(figsize=(4, 3))
sns.barplot(data=cog_counts, y="Function", x="Count", color="skyblue")
plt.xlabel("Number of Proteins")
plt.ylabel("COG Functional Category")
#plt.title("COG Categories in Widespread Plasmid Proteins")
plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_plasmids_cog_category_plot_filtered.png", dpi=300)
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_plasmids_cog_category_plot_filtered.svg")

plt.show()
