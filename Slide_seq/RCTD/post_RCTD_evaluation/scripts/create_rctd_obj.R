# bash
# apptainer shell \
#   --bind /scratch/mfafouti/miniforge3/envs/spacexr_env:/opt/spacexr_env \
#   --bind /scratch/mfafouti/miniforge3/envs/spacexr_env/lib/R/library:/opt/spacexr_env/lib/R/library \
#   --bind /scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run:/workspace \
#   --bind /scratch/mfafouti/Mommybrain/Slide_seq/RCTD/post_RCTD_evaluation/scripts:/workspace/scripts \
#   --bind /scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/B03:/workspace/out_RCTD_all/B03 \
#   --bind /scratch/mfafouti/Mommybrain/Slide_seq/RCTD/post_RCTD_evaluation/RCTD_object_QC:/workspace/RCTD_object_QC \
#   /scratch/mfafouti/Mommybrain/seurat_v5.sif 


library(spacexr, lib.loc = "/opt/spacexr_env/lib/R/library")
library(dplyr)
library(tidyr)
library(purrr)

# PATHS
# sample_output = "/workspace/out_RCTD_all/B03"
# rds_files <- list.files(sample_output, pattern = "B03_subset_.*_RCTD\\.rds$", full.names = TRUE)
rds_file <- "/workspace/out_RCTD_all/B03/B03_subset_1_RCTD.rds"

# --- Read RDS ---
rctd_obj <- readRDS(rds_file)

# --- Extract the first subset ---
subset_vars <- rctd_obj@spatialRNA

names(subset_vars)
str(subset_vars)





# head(subset_vars, 5)

# # Inspect results slot (usually contains per-cell info)
# names(rctd_file@results)
# str(rctd_file@results)

# # Optional: see number of cells and predicted weights
# dim(rctd_file@results$weights)   # cells x reference cell types
# head(rctd_file@results$weights)

# # Check assigned cell types
# head(rctd_file@results$cell_type)

# # Check doublets if present
# head(rctd_file@results$doublet_score)



# internal_vars_df <- imap_dfr(rctd_list, function(obj, idx) {
#   iv <- obj@internal_vars
  
  # Flatten each element of the list
#   map2_dfr(names(iv), iv, function(varname, vec) {
#     tibble(
#       variable = varname,
#       value = as.character(vec),   # convert everything to character
#       file = basename(rds_files[idx])
#     )
#   })
# })

# internal_vars_wide <- imap_dfr(rctd_list, function(obj, idx) {
#   iv <- obj@internal_vars
#   df <- map(iv, function(x) {
#     if(length(x) > 1) paste(x, collapse = ",") else as.character(x)
#   }) %>% as_tibble()
#   df$file <- basename(rds_files[idx])
#   df
# })

# # Save wide CSV
# write.csv(internal_vars_wide, "/workspace/RCTD_object_QC/internal_vars_summary_wide.csv", row.names = FALSE)

# Save to CSV
#write.csv(internal_vars_df, "/workspace/RCTD_object_QC/1subset_internal_vars_summary.csv", row.names = FALSE)