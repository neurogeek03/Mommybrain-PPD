# =====================================================================
# Title:        Convert to H5Seurat
# Description:  Using a .rds file as input, a .h5seurat file is created.
# Author:       Maria Eleni Fafouti
# Date:         2025-04-25
# =====================================================================

# bash
# conda activate seurat_env

# ---- 1. SETUP ----
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)

# Define input and output file paths
project_path <- "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags"
in_path <- file.path(project_path, 'Integration', 'Integration_CB_out') 
doublet_folder <- file.path(project_path, 'Doublet_detection')

# ---- 2. DATA LOADING & PRE-PROCESSING ----
# Read in your Seurat object
seu <- readRDS(file.path(in_path, "new_integrated.rds"))
message("Seurat object loaded!")

# Checking the contents of the object
head(seu@meta.data)
table(seu$orig.ident)

# Joining Layers 
seu[["RNA"]] <- JoinLayers(seu[["RNA"]]) 
message("Layers Joined!")

# Switching to the RNA assay (raw data) - not the integrated one - because DoubletFinder requires unintegrated data
DefaultAssay(seu) <- "RNA"
message("Assay switched!")

# ---- 3. DATA LOADING & PRE-PROCESSING ----
# Converting to SingleCellExperiment
sce <- as.SingleCellExperiment(seu)

# ---- 4. DATA LOADING & PRE-PROCESSING ----
saveRDS(sce, file.path(doublet_folder, "sce_object.rds"))

