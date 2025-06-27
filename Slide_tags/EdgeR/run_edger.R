# ==========================================================================
# Title:        Run DE using EdgeR
# Description:  Running EdgeR LRT model and saving the results as csv files 
# Author:       Maria Eleni Fafouti
# Date:         2025-06-18
# ==========================================================================

# bash 
# cd to mommybrain folder on scratch

# apptainer shell \
#   --bind /gpfs/fs0/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/EdgeR:/workspace \
#   --bind /gpfs/fs0/scratch/s/shreejoy/mfafouti/miniforge3/envs/rctd_env:/opt/rctd_env \
#   --bind /gpfs/fs0/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Reference_data/new:/project \
#   /gpfs/fs0/scratch/s/shreejoy/mfafouti/Mommybrain/edger.sif

# ========= LIBRARIES =========
library(edgeR)
library(stringr, lib.loc = "/opt/rctd_env/lib/R/library")
library(vctrs, lib.loc = "/opt/rctd_env/lib/R/library")

# ========= CONTAINER PATHS =========
input_dir <- "/workspace/pseudobulk_output_25k_genes"
output_dir <- "/workspace/edger_out"

# List all counts files
counts_files <- list.files(input_dir, pattern = "_counts.tsv$", full.names = TRUE)

for (counts_file in counts_files) {
  celltype <- str_replace(basename(counts_file), "_counts.tsv", "")
  message("🔬 Processing cell type: ", celltype)
  
  # Load data
  counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1)
  meta_file <- file.path(input_dir, "metadata_universal.tsv")
  metadata <- read.table(meta_file, header = TRUE, sep = "\t", row.names = 1)
  
  print(colnames(counts))
  print(rownames(metadata))

  # Check order
  stopifnot(all(colnames(counts) == rownames(metadata)))

  # edgeR analysis
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  
  metadata$treatment <- factor(metadata$treatment)
  design <- model.matrix(~ treatment, data = metadata)
  
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = 2)  # Assumes treatmentB is second level
  
  # Save results
  result_file <- file.path(output_dir, paste0(celltype, "_dge_results.tsv"))
  write.table(topTags(lrt, n = Inf)$table, result_file, sep = "\t", quote = FALSE)
}

message("✅ edgeR analysis complete for all cell types.")