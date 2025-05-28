# =====================================================================
# Title:        Post Integration Analysis 
# Description:  Performing dimensionality reduction computations on integrated data.
# Author:       Maria Eleni Fafouti
# Date:         2025-04-20
# =====================================================================

# ---- 1. SETUP ----
library(Seurat)
library(dplyr)
library(ggplot2)

rds_path = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Integration_CB_out/new_integrated.rds"
output_dir = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_integration_out"
integrated <- readRDS(rds_path)

# Preview data
integrated
head(integrated@meta.data)
unique(integrated@meta.data$orig.ident)
unique(integrated@meta.data$treatment)

# ---- 2. DIMENSIONALITY REDUCTION ----
integrated <- integrated %>%
    #JoinLayers(integrated) %>%
    ScaleData() %>%
    RunPCA() 

# saving the elbow plot to diagnose PCA
elbow_cb_int <- ElbowPlot(integrated, ndims = 50)
elbow_path = file.path(output_dir, "elbow_cb_int.png")
ggsave(filename = elbow_path, plot = elbow_cb_int, width = 8, height = 6, dpi = 300)

message("Elbow plot saved at: ", elbow_path)

# finding neighbors, clustering, UMAP
integrated <- integrated %>%
    FindNeighbors(dims = 1:30) %>% #Using the first 30 PCs
    FindClusters() %>%
    RunUMAP(dims = 1:30)
    # RunUMAP(graph = "integrated_snn") #this is already what happens by defalut

message("Computations completed!")

# ---- 3. PLOTTING UMAPS ----
# Generate and save the UMAP plot
# Plotting the umap & coloring based on sample treatment
groupings <- c("treatment", "seurat_clusters", "orig.ident")
titles <- c("Colored by Treatment", "Colored by Seurat Clusters", "Colored by Sample")
file_names <- c("umap_treatment.png", "umap_clusters.png", "umap_samples.png")

# Loop through each grouping and create the UMAP plot
for (i in 1:length(groupings)) {
  umap_plot <- DimPlot(integrated, reduction = "umap", group.by = groupings[i]) +
    ggtitle(paste("UMAP of Integrated Seurat Object -", titles[i]))
  
  output_file <- paste(output_dir, file_names[i], sep = "")
  
  ggsave(output_file, plot = umap_plot, width = 8, height = 6, dpi = 300)
  
  message(paste("UMAP plot saved as", file_names[i]))
}

message("UMAP plots colored by", groupings,"were plotted and saved!", "\n")
