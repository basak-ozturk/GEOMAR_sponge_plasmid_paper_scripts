# ===============================
# Top 10 KEGG Pathway Bubble Plot
# ===============================

library(tidyverse)
library(scales)
library(svglite)
library(stringr)

# --- Load data ---
data <- read.table(
  "C:/Users/hayat/Downloads/R_files/data/KEGG_Pathway_per_plasmid_with_names.txt",
  header = FALSE,
  sep = "\t",
  col.names = c("Plasmid", "Pathways")
)

# --- Process data ---
data_expanded <- data |>
  separate_rows(Pathways, sep = ",") |>
  count(Pathways, name = "Count") |>
  filter(Count >= 100) |>
  arrange(desc(Count)) |>
  slice(1:10)  # Keep only top 10 pathways

# --- Wrap long pathway names to reduce vertical space ---
data_expanded$Pathways <- str_wrap(data_expanded$Pathways, width = 25)

# --- Compute dynamic x-axis breaks ---
x_range <- range(log10(data_expanded$Count))
x_breaks <- pretty(x_range, n = 6)  # ~5–6 ticks

# --- Define plot ---
p <- ggplot(
  data_expanded,
  aes(x = log10(Count), y = reorder(Pathways, Count), size = Count)
) +
  geom_point(alpha = 0.8, color = "steelblue", stroke = 0.4) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = function(x) round(x, 1)
  ) +
  scale_size_continuous(trans = "log10", range = c(2, 8)) +
  scale_y_discrete(expand = expansion(mult = 0.05)) +
  labs(
    x = "log10(Number of Plasmids)",
    y = "Pathway",
    title = "Top 10 Pathway Enrichment (≥100 Plasmids)",
    size = "Plasmid Count"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    axis.text.y = element_text(size = 8, face = "bold", hjust = 1),
    axis.text.x = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    title = element_text(size = 8),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  coord_cartesian(clip = "off")  # prevent text cut-off

# --- Display plot for preview ---
print(p)

# --- Estimate figure height (0.35 in per pathway) ---
n_pathways <- nrow(data_expanded)
height_in <- max(4, n_pathways * 0.35)  # minimum 4 in

# --- Save as SVG (publication-ready) ---
ggsave(
  "C:/Users/hayat/Downloads/R_files/graphs/KEGG_Pathway_BubblePlot_with_names_top10.svg",
  plot = p,
  device = svglite,
  width = 4,        # 100 mm
  height = height_in,
  units = "in",
  bg = "white"
)
