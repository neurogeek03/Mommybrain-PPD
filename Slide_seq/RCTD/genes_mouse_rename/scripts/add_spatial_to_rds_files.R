library(Seurat)
library(dplyr)

# Paths
rds_dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/UPDATED_mouse_query_rds"
csv_dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/UPDATED_collapsed_adata_with_mouse_orthologs/spatial_info_csvs"
output_dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/UPDATED_mouse_query_rds_files_w_spatial"
dir.create(output_dir, showWarnings = FALSE)

# List all .rds files
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

for (rds_file in rds_files) {
  # Load the RDS file (Seurat object)
  seu <- readRDS(rds_file)
  
  # Get the sample name (before first '_')
  sample_name <- strsplit(basename(rds_file), "\\.")[[1]][1]

  print(sample_name)
  
  # Find matching CSV file
  csv_file <- file.path(csv_dir, paste0(sample_name, "_spatial_coords.csv"))
  
  if (!file.exists(csv_file)) {
    message("No spatial CSV found for ", sample_name, ", skipping.")
    next
  }
  
  # Read spatial coordinates CSV
  coords <- read.csv(csv_file, row.names = 1)
  
  # Make sure cell names match
  common_cells <- intersect(rownames(seu@meta.data), rownames(coords))
  if (length(common_cells) == 0) {
    warning("No matching cells between Seurat object and CSV for ", sample_name)
    next
  }
  
  # Add to metadata
  seu@meta.data[common_cells, "spatial_x"] <- coords[common_cells, "x"]
  seu@meta.data[common_cells, "spatial_y"] <- coords[common_cells, "y"]
  
  # Save updated object
  saveRDS(seu, file = file.path(output_dir, paste0(sample_name, "_with_spatial.rds")))
  
  message("Updated ", sample_name, " with spatial metadata.")
}
