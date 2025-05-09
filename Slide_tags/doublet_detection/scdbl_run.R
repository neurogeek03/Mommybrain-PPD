# =====================================================================
# Title:        Run scDblFinder on H5AD Files and Export Doublet Labels
# Description:  Iterates over h5ad files, runs scDblFinder, and saves
#               doublet classifications and scores to CSV.
# Author:       Maria Eleni Fafouti
# Date:         2025-04-25
# =====================================================================

# ---- SETUP ----

# Ensure you have activated the correct conda environment before running:
# conda activate scdbl_env

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(zellkonverter)
  library(dplyr)
  library(readr)
})

sample_list <- c("BC28", "BC3", "BC9", "BC15", "BC14", "BC13")

# ---- DEFINE PATHS ----
project_path <- "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Doublet_detection/scanpy"
in_path <- file.path(project_path, "tmp_dir")

# ---- PROCESS EACH SAMPLE ----
for (sample in sample_list) {
  message(paste("Processing sample:", sample))

  sample_file <- paste0(sample, "_tmp.h5ad")
  out_csv_name <- paste0(sample, "_scdblfinder.csv")
  in_file <- file.path(in_path, sample_file)
  csv_file <- file.path(in_path, out_csv_name)

  # ---- LOAD DATA ----
  sce <- readH5AD(in_file)

  if (!"counts" %in% assayNames(sce)) {
    assay(sce, "counts") <- assay(sce, "X")
    assays(sce) <- assays(sce)["counts"]
  }

  # ---- RUN SCDBLFINDER ----
  sce <- scDblFinder(sce, BPPARAM = SerialParam())

  # ---- SAVE RESULTS TO CSV ----
  label_df <- data.frame(
    barcode = colnames(sce),
    scDblFinder.class = sce$scDblFinder.class,
    scDblFinder.score = sce$scDblFinder.score
  )

  write_csv(label_df, csv_file)
}

