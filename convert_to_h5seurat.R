# =====================================================================
# Title:        Convert to H5Seurat
# Description:  Using a .rds file as input, a .h5seurat file is created.
# Author:       Maria Eleni Fafouti
# Date:         2025-03-27
# =====================================================================

# ---- 1. SETUP ----
# Load required packages
library(Seurat)
library(SeuratDisk)

# Define input and output file paths
input_rds <- "/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/OUTPUTS2/results/OUTPUT/B21/B21_seurat.rds"
output_h5seurat <- "/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/background_removal/B21/B21_seurat.h5seurat"

# ---- 2. DATA LOADING ----
# Read the RDS file 
seurat_obj <- readRDS(input_rds)
cat("✅ Object loaded from:", input_rds, "\n")

# ---- 3. OUTPUT ----
# Save as H5Seurat
SaveH5Seurat(seurat_obj, filename = output_h5seurat, overwrite = TRUE)
cat("✅ H5Seurat file saved to:", output_h5seurat, "\n")

# ---- END OF SCRIPT ----