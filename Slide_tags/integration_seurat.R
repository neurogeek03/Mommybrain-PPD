# =====================================================================
# Title:        Seurat Integration Analysis
# Description:  Integrate different Seurat Objects and create a UMAP plot
# Author:       Maria Eleni Fafouti
# Date:         2025-04-02
# =====================================================================

# This script is submitted as a job, accoring to the sbatch_integration.sh file 

# ---- 1. SETUP ----
# Load required packages
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Define the output base directory
output_dir <- "/home/mfafouti/scratch/Mommybrain_marlen/Slide_tags/Single_cell_analysis/Integration_out/"

# ---- 2. DATA LOADING ----
# Automatically getting the rds paths 
base_dir <- "/home/mfafouti/nearline/def-shreejoy/RatAllCTOuts"
subdirs <- list.dirs(base_dir, recursive = FALSE)
# Filter: Keep only subdirectories where the numeric part before "_" ends in "1"
subdirs <- subdirs[grepl("/[0-9]*1_", subdirs)]
seurat_list <- list()

# Loop through each subdirectory
for (subdir in subdirs) {
  # Extract the sample name (last part after '_')
  sample_name <- sub(".*_", "", basename(subdir))

  # Construct the expected file path
  rds_path <- file.path(subdir, paste0("trekker_", sample_name), "output", paste0(sample_name, "_ConfPositioned_seurat_spatial.rds"))

  # Check if the file exists before reading
  if (file.exists(rds_path)) {
    seurat_list[[sample_name]] <- readRDS(rds_path)
    message("Loaded: ", sample_name)
  } else {
    message("File not found: ", rds_path)
  }
}

cat("✅ Objects loaded from base directory:", base_dir, "\n")

# ---- 3. INTEGRATION OF SAMPLES ----
# Identify integration features
features <- SelectIntegrationFeatures(object.list = seurat_list)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)

# Integrate datasets
integrated <- IntegrateData(anchorset = anchors)

# Continue with standard Seurat workflow
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

# ---- 4. SAVING INTEGRATED OBJECT ----
# Define output file paths
rdata_path <- file.path(output_dir, "integrated_seurat.RData")
umap_path <- file.path(output_dir, "integrated_umap.png")

# Save the integrated Seurat object
save(integrated, file = rdata_path)
message("Integrated Seurat object saved at: ", rdata_path)

# ---- 5. PLOTTING UMAP ----
# Generate and save the UMAP plot
umap_plot <- DimPlot(integrated, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP of Integrated Seurat Object")

ggsave(umap_path, plot = umap_plot, width = 8, height = 6, dpi = 300)
message("UMAP plot saved at: ", umap_path)

# # ---- END OF SCRIPT ----