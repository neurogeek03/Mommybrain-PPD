#!/usr/bin/env Rscript

# Usage: Rscript export_rctd_weights.R <path_to_rctd_object.rds> <output_directory>

library(spacexr, lib.loc = "/opt/miniforge3/envs/rctd_env/lib/R/library")
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript export_rctd_weights.R <path_to_rctd_object.rds> <output_directory>")
}

rctd_path <- args[1]
out_dir   <- args[2]

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("Loading RCTD object from:", rctd_path, "
")
rctd <- readRDS(rctd_path)

# 1. Export summary results (first_type, second_type, weights, etc.)
cat("Exporting summary results...
")
if (!is.null(rctd@results$results_df)) {
  write.csv(rctd@results$results_df, 
            file = file.path(out_dir, "rctd_summary_results.csv"), 
            row.names = TRUE)
}

# 2. Export full weight matrix (all cell types for all pixels)
cat("Exporting full weight matrix...
")
if (!is.null(rctd@results$weights)) {
  # weights is typically a sparse matrix
  weights_mat <- as.matrix(rctd@results$weights)
  write.csv(weights_mat, 
            file = file.path(out_dir, "rctd_full_weights.csv"), 
            row.names = TRUE)
  
  # 3. Export normalized weights (proportions summing to 1 per pixel)
  cat("Exporting normalized weights...
")
  norm_weights <- sweep(weights_mat, 1, rowSums(weights_mat), "/")
  # Handle potential division by zero for pixels with no detected signal
  norm_weights[is.na(norm_weights)] <- 0
  write.csv(norm_weights, 
            file = file.path(out_dir, "rctd_normalized_weights.csv"), 
            row.names = TRUE)
}

cat("Done. Results saved to:", out_dir, "
")
