# =====================================================================
# Title:        Convert to H5Seurat
# Description:  Using a .rds file as input, a .h5seurat file is created.
# Author:       Maria Eleni Fafouti
# Date:         2025-04-25
# =====================================================================

# bash
# conda activate scdbl_env

# ---- 1. SETUP ----

# Library Loading
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)

# Define input and output file paths
project_path <- "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags"
in_path <- file.path(project_path, 'Integration', 'Integration_CB_out') 
doublet_folder <- file.path(project_path, 'Doublet_detection')

# ---- 2. DATA LOADING ----
# Loading the SingleCellExperiment obj
sce <- readRDS(file.path(doublet_folder, 'sce_object.rds'))

# Preview 
sce

# ---- 3. RUNNING SCDBL FINDER ----
sce <- scDblFinder(sce, samples="orig.ident", BPPARAM = SerialParam())

# Verifying 
table(sce$scDblFinder.class)

# ---- 4. SAVING DATA AS SINGLE CELL EXPERIMENT OBJECT ----
saveRDS(sce, file.path(doublet_folder, "scdblfinder_result1.rds"))

