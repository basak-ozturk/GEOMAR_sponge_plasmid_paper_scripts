from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import pandas as pd

# Define input files
file1 = "C:/Users/hayat/Downloads/R_files/data/plasmids_with_plasmid_hallmark.txt"
file2 = "C:/Users/hayat/Downloads/R_files/data/oriT_alloriT_blast_results_names.txt"
file3 = "C:/Users/hayat/Downloads/R_files/data/unique_CONJscan_hit_ids.txt"

# Read content into sets
with open(file1, "r") as f:
    set1 = set(line.strip() for line in f if line.strip())

with open(file2, "r") as f:
    set2 = set(line.strip() for line in f if line.strip())

with open(file3, "r") as f:
    set3 = set(line.strip() for line in f if line.strip())

# --- Compute all Venn subsets ---
output = {
    "only_set1": sorted(set1 - set2 - set3),
    "only_set2": sorted(set2 - set1 - set3),
    "only_set3": sorted(set3 - set1 - set2),
    "set1_set2": sorted((set1 & set2) - set3),
    "set1_set3": sorted((set1 & set3) - set2),
    "set2_set3": sorted((set2 & set3) - set1),
    "set1_set2_set3": sorted(set1 & set2 & set3)
}

# Convert to a long-table CSV
rows = []
for category, items in output.items():
    for item in items:
        rows.append({"category": category, "id": item})

df = pd.DataFrame(rows)
df.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/venn3_output.csv",
    index=False
)

print("CSV saved with all Venn memberships.")


# --------- Plot the Venn Diagram ---------
plt.figure(figsize=(8, 8))
venn = venn3(
    [set1, set2, set3],
    set_labels=("Plasmid Hallmark", r"$\it{oriT}$", "CONJscan"),
    set_colors=("#E63946", "#457B9D", "#2A9D8F"),
    alpha=1
)

plt.show()
plt.close()
