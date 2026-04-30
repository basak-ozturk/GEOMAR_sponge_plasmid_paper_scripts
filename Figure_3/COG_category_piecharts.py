import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
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

TOTAL = 113743  # your full input protein count

df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations", sep="\t", comment="#")

# Normalize missing values
for col in ["COG_category", "PFAMs", "KEGG_Pathway"]:
    df[col] = df[col].replace("-", pd.NA)

# Counts
n_output = len(df)  # proteins that got ANY annotation line from eggnog
n_no_output = TOTAL - n_output

n_cog = df["COG_category"].notna().sum()
n_cog_missing = TOTAL - n_cog

def has_informative_cog(x):
    if pd.isna(x): return False
    return any(letter not in ["R", "S"] for letter in str(x))

n_cog_inform = df["COG_category"].apply(has_informative_cog).sum()
n_cog_uninform = TOTAL - n_cog_inform

n_pfam = df["PFAMs"].notna().sum()
n_pfam_missing = TOTAL - n_pfam

n_kegg = df["KEGG_Pathway"].notna().sum()
n_kegg_missing = TOTAL - n_kegg

# Build pie chart data
categories = [
    ("EggNOG output", n_output, n_no_output),
    ("COG category", n_cog, n_cog_missing),
    ("Informative COG (not R/S)", n_cog_inform, n_cog_uninform),
    ("PFAM", n_pfam, n_pfam_missing),
    ("KEGG pathway", n_kegg, n_kegg_missing),
]

# Plot
fig, axes = plt.subplots(2, 3, figsize=(7, 4.5))
axes = axes.flatten()

for ax, (label, assigned, missing) in zip(axes, categories):
    ax.pie(
        [assigned, missing],
        labels=[f"Assigned ({assigned})", f"Not assigned ({missing})"],
        autopct='%1.1f%%',
        colors=["steelblue", "lightgray"]
    )
    ax.set_title(label)

# remove last empty subplot if present
axes[-1].axis("off")

plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/eggnog_annotation_pies.png", dpi=300)
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/eggnog_annotation_pies.svg")

plt.show()

labels = [
    "EggNOG output",
    "COG",
    "Informative COG",
    "PFAM",
    "KEGG"
]

assigned = [
    n_output,
    n_cog,
    n_cog_inform,
    n_pfam,
    n_kegg
]

missing = [
    n_no_output,
    n_cog_missing,
    n_cog_uninform,
    n_pfam_missing,
    n_kegg_missing
]

x = np.arange(len(labels))

plt.figure(figsize=(5,3))
plt.bar(x, assigned, label="Assigned", color="steelblue")
plt.bar(x, missing, bottom=assigned, label="Not assigned", color="lightgray")

plt.xticks(x, labels, rotation=45, ha='right')
plt.ylabel("Proteins")
plt.title("Annotation completeness per category")
plt.legend()
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/annotation_stacked_bar.png", dpi=300)
plt.show()

