#!/usr/bin/env Rscript

# -------------------------------
# Marker discovery + dot plot script
# -------------------------------

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)


# -------------------------------
# Load your Seurat object
# -------------------------------
# Example: load from RDS
pbmc <- readRDS("/workspace/data/slide_tags_reference_merged.rds")

print("object loaded!")

pbmc <- FindNeighbors(pbmc, dims = 1:30)  # adjust dims as appropriate
pbmc <- FindClusters(pbmc, resolution = 0.5)
table(Idents(pbmc))
# -------------------------------
# Find markers
# -------------------------------
all_markers <- FindAllMarkers(object = pbmc) 
print("The column names of all markers are:")
colnames(all_markers)

all_markers <- all_markers%>%
  Add_Pct_Diff() %>%
  filter(pct_diff > 0.6)

print("markers found!")
# -------------------------------
# Extract top markers
# -------------------------------
top_markers <- Extract_Top_Markers(
  marker_dataframe = all_markers,
  num_features = 7,
  named_vector = FALSE,
  make_unique = TRUE
)
print("top markers extracted!")
# -------------------------------
# Plot markers
# -------------------------------
p <- Clustered_DotPlot(
  seurat_object = pbmc,
  features = top_markers
)
print("plot created!")

# Save the plot
ggsave("clustered_dotplot.png", p, width = 10, height = 8, dpi = 300)

print("plot saved!")
# -------------------------------
# Print session info
# -------------------------------
sessionInfo()
