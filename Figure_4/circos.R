library(circlize)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(Cairo)   # <-- important for editable SVG output

# === Load data ===
nodes <- read.csv("C:/Users/hayat/Downloads/R_files/data/nodes_with_host_groups.csv")
host_map <- read.csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_hosts.csv")
mcl_df <- read.csv("C:/Users/hayat/Downloads/R_files/data/cytoscape_edges_with_mcl_clusters.csv")

# Extract unique node-cluster relationships
mcl_nodes <- mcl_df |>
  select(name, X__mclCluster) |>
  rename(ID = name, Cluster = X__mclCluster) |>
  distinct()

# Merge all data
merged_df <- nodes |>
  inner_join(host_map, by = "ID") |>
  left_join(mcl_nodes, by = "ID")

# Count plasmids per Host ↔ Host_Group
chord_data <- merged_df |>
  group_by(Host, Host_Group) |>
  summarise(Count = n(), .groups = "drop")

# Convert to wide matrix
chord_matrix <- pivot_wider(
  chord_data,
  names_from = Host_Group,
  values_from = Count,
  values_fill = 0
)

colnames(chord_matrix) <- gsub("[–—−]", "-", colnames(chord_matrix))
desired_order <- c("1", "2-3", "4-7", "8-14", "15+")
chord_matrix <- chord_matrix[, desired_order]
all_hosts <- sort(unique(merged_df$Host))
rownames(chord_matrix) <- all_hosts
chord_matrix <- as.matrix(chord_matrix)

# === Colors ===
set.seed(42)
host_colors <- setNames(viridis(length(all_hosts), option = "C"), all_hosts)
all_host_groups <- sort(unique(nodes$Host_Group))
set.seed(24)
host_group_colors <- setNames(brewer.pal(length(all_host_groups), "Set3"), all_host_groups)
grid.col <- c(host_colors, host_group_colors)

# === Editable SVG output ===
CairoSVG(
  file = "C:/Users/hayat/Downloads/R_files/graphs/circos_plot_host_count_editable.svg",
  width = 3.15,   # ≈ 80 mm
  height = 3.15,
  pointsize = 8
)

par(family = "Arial")  # use a clean font for Inkscape text editing

chordDiagram(
  chord_matrix,
  grid.col = grid.col,
  transparency = 0.2,
  annotationTrack = c("grid", "name"),  # keep host & group labels
  preAllocateTracks = 1
)

dev.off()
