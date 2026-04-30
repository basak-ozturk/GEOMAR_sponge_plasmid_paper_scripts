import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import os
#import matplotlib.image as mpimg
#import numpy as np
#import matplotlib.colors as mcolors
#import matplotlib.cm as cm
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

# === File Paths ===
data_folder = "C:/Users/hayat/Downloads/R_files/data/"
output_folder = "C:/Users/hayat/Downloads/R_files/graphs/"
os.makedirs(output_folder, exist_ok=True)

location_file = os.path.join(data_folder, "top_widespread_plasmid_metadata.tsv")
density_file = os.path.join(data_folder, "plasmid_density_summary.csv")
cluster_file = os.path.join(data_folder, "cytoscape_edges_with_mcl_clusters.csv")  # <- update if needed

# === 1. Load Data ===
location_df = pd.read_csv(location_file, sep="\t")
density_df = pd.read_csv(density_file)
cluster_df = pd.read_csv(cluster_file)

# === 2. Prepare cluster labels ===
cluster_df = cluster_df.rename(columns={"name": "Plasmid", "__mclCluster": "Cluster"})
cluster_df["Cluster"] = cluster_df["Cluster"].apply(lambda x: f"Cluster {int(x)}" if pd.notna(x) else "Singleton")

# Merge cluster info into location data
location_df = location_df.merge(cluster_df[["Plasmid", "Cluster"]], on="Plasmid", how="left")
location_df["Cluster"] = location_df["Cluster"].fillna("Singleton")

# === 3. Assign location IDs by rounded coords ===
bin_size = 2  # degrees, change to higher you want coarser

location_df["Lat_round"] = (location_df["Latitude"] // bin_size) * bin_size
location_df["Lon_round"] = (location_df["Longitude"] // bin_size) * bin_size

location_df["Location_ID"] = location_df.groupby(["Lat_round", "Lon_round"]).ngroup() + 1

# Get unique locations with coordinates and IDs for map plotting
unique_locations = location_df.drop_duplicates(subset=["Lat_round", "Lon_round", "Location_ID"])[
    ["Lat_round", "Lon_round", "Location_ID"]].reset_index(drop=True)

# === 4. Plot Density + Numbered Locations Map (Plotly) ===

# Density background trace
density_trace = go.Scattergeo(
    lat=density_df["Latitude_round"],
    lon=density_df["Longitude_round"],
    mode="markers",
    marker=dict(
        size=density_df["Plasmid_Count"],
        color=density_df["Plasmid_Count"],
        colorscale='Greys',
        cmin=density_df["Plasmid_Count"].min(),
        cmax=density_df["Plasmid_Count"].max(),
        colorbar=dict(title="Plasmid Count"),
        line=dict(width=1, color='red'),
        sizemode='area',
        sizeref=2. * density_df["Plasmid_Count"].max() / 25 ** 2,
        sizemin=2
    ),
    name="Plasmid Density",
    hovertemplate="Plasmid Count: %{marker.size}<extra></extra>"
)

# Location markers with numbers
location_markers = go.Scattergeo(
    lat=unique_locations["Lat_round"],
    lon=unique_locations["Lon_round"],
    mode="markers+text",
    marker=dict(
        size=6,
        color="blue",
        symbol="circle",
        line=dict(width=1, color="black")
    ),
    text=unique_locations["Location_ID"].astype(str),
    textposition="top center",
    name="Sampling Locations",
    hoverinfo="text",
    hovertext=["Location ID: " + str(i) for i in unique_locations["Location_ID"]],
)

# Build figure
# Set up Plotly figure
fig = go.Figure(data=[density_trace, location_markers])
fig.update_layout(
    title="Plasmid Density Map with Numbered Locations",
    geo=dict(
        showland=True,
        landcolor="rgb(243, 243, 243)",
        projection_type="natural earth",
        showcountries=True,
    ),
    margin=dict(r=0, t=40, l=0, b=0),
    legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
)

# Save as SVG with width ~80 mm
map_svg_path = os.path.join(output_folder, "plasmid_density_numbered_map.svg")
#fig.write_image(map_svg_path, width=800, height=600)  # width px ~ 80 mm at 300 dpi
print(f"Map saved as SVG: {map_svg_path}")


# === 5. Create Pie Charts per Location (Matplotlib) ===

cluster_counts = location_df.groupby(["Location_ID", "Cluster"]).size().reset_index(name="Count")
unique_clusters = sorted(location_df["Cluster"].unique())

color_palette = px.colors.qualitative.Dark24 + px.colors.qualitative.Pastel + px.colors.qualitative.Set3
color_map = {cluster: color_palette[i % len(color_palette)] for i, cluster in enumerate(unique_clusters)}

# === Determine plasmid counts per location ===
location_counts = (
    location_df
    .groupby("Location_ID")
    .size()
    .reset_index(name="Plasmid_Count")
)

# === Plot pie charts for multi-plasmid locations ===
filtered_locations = location_counts[location_counts["Plasmid_Count"] >= 1]
# Get all unique clusters
unique_clusters = sorted(cluster_counts["Cluster"].unique())

extended_colors = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
    "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7",
    "#dbdb8d", "#9edae5", "#393b79", "#637939", "#8c6d31"
]
color_map = {cluster: extended_colors[i % len(extended_colors)] for i, cluster in enumerate(unique_clusters)}

# Adjust figure size for 80 mm width
width_inch = 90 / 25.4  # 80 mm → inches
cols = 5
rows = (len(filtered_locations) + cols - 1) // cols
fig_pies, axes = plt.subplots(rows, cols, figsize=(width_inch, rows * width_inch / cols), subplot_kw=dict(aspect="equal"))
axes = axes.flatten()

for ax, (_, row) in zip(axes, filtered_locations.iterrows()):
    loc_id = row["Location_ID"]
    plasmid_count = row["Plasmid_Count"]
    data = cluster_counts[cluster_counts["Location_ID"] == loc_id]
    if data.empty:
        ax.axis('off')
        continue
    sizes = data["Count"]
    labels = data["Cluster"]
    colors = [color_map[c] for c in labels]
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=None,
        colors=colors,
        autopct='%1.1f%%',
        startangle=140,
        textprops={'fontsize': 4}
    )
    ax.set_title(f"Location {loc_id} ({plasmid_count} plasmids)", fontsize=8)

# Hide unused axes
for j in range(len(filtered_locations), len(axes)):
    axes[j].axis('off')

# Add legend below the figure in 3 columns
fig_pies.legend(
    unique_clusters,
    title="Clusters",
    loc="lower center",
    bbox_to_anchor=(0.5, -0.05),  # move below figure
    fontsize=4,
    ncol=3,  # number of columns in legend
    frameon=False  # optional, for cleaner look
)

plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave space at bottom for legend

# Save as SVG
pie_svg_path = os.path.join(output_folder, "cluster_composition_piecharts_all.svg")
fig_pies.savefig(pie_svg_path, format="svg", dpi=300, bbox_inches="tight")
plt.close(fig_pies)
print(f"Pie charts saved as SVG: {pie_svg_path}")


# Prepare proportion table
composition_df = (
    cluster_counts
    .pivot(index="Location_ID", columns="Cluster", values="Count")
    .fillna(0)
)

# Normalize to proportions
composition_df = composition_df.div(composition_df.sum(axis=1), axis=0)

# Make figure ~60–75% of current size
fig, ax = plt.subplots(figsize=(3.5, 3.6))  # was (6, 4)

# Plot stacked bars
composition_df.plot(
    kind="barh",
    stacked=True,
    color=[color_map[c] for c in composition_df.columns],
    ax=ax
)

# Labels and legend
ax.set_xlabel("Proportion of ANI Clusters", fontsize=8)
ax.set_ylabel("Location", fontsize=8)
ax.tick_params(axis='both', which='major', labelsize=7)

ax.legend(
    #title="Cluster",
    bbox_to_anchor=(0.5, -0.25),
    loc="upper center",
    fontsize=7,
    title_fontsize=7,
    ncol=3,  # spread across several columns
    frameon=False
)


plt.tight_layout()
plt.savefig(
    os.path.join(output_folder, "cluster_composition_stackedbars.svg"),
    bbox_inches="tight",
    dpi=300
)
plt.show()



output_csv_path = "C:/Users/hayat/Downloads/R_files/data/location_df_saved.csv"
location_df.to_csv(output_csv_path, index=False)



