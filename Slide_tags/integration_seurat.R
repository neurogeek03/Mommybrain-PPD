# =====================================================================
# Title:        Seurat Integration Analysis
# Description:  Integrate different Seurat Objects and create a UMAP plot
# Author:       Maria Eleni Fafouti
# Date:         2025-04-14
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
output_dir <- "/home/mfafouti/scratch/Mommybrain_marlen/Slide_tags/Single_cell_analysis/Integration/Integration_out"

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

seurat_list[["BC28"]] <- readRDS("/home/mfafouti/scratch/Mommybrain_marlen/Slide_tags/trekker_BC28/output/BC28_ConfPositioned_seurat_spatial.rds")  
print("✅ Object BC28 loaded!")

# ---- 3. INTEGRATION OF SAMPLES ----

seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  return(obj)
})
print("✅ Normalization!")

# Identify integration features
features <- SelectIntegrationFeatures(object.list = seurat_list)
print("✅ Selection of Integration features!")

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
print("✅ Finding anchors!")

# Integrate datasets
integrated <- IntegrateData(anchorset = anchors)
print("✅ Integration done!")

# Setting the integrated assay as the default object 
DefaultAssay(integrated) <- "integrated"

# Adding a metadata column to show what treatment each sample received
integrated@meta.data$treatment <- ifelse(
  integrated@meta.data$orig.ident %in% c("BC3", "BC9", "BC15"),
  "CORT",
  "OIL"
)

# ---- 4. SAVING THE INTEGRATED OBJECT ----

# As an .rdata file
rds_path <- file.path(output_dir, "new_integrated.rds")
saveRDS(integrated, file = rds_path)

# As a .rds file
rdata_path <- file.path(output_dir, "new_integrated.RData")
save(integrated, file = rdata_path) #TODO save this as a integration script and push to repo
print("✅ Object saved!")

# ---- END OF SCRIPT ----

