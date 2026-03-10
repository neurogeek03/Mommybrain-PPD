# ==========================================================================
# Title:        Run DE using EdgeR
# Description:  Running EdgeR LRT model and saving the results as csv files 
# Author:       Maria Eleni Fafouti
# Date:         2025-06-18
# ==========================================================================

# bash 
# cd to mommybrain folder on scratch

# apptainer shell \
#   --bind /scratch/mfafouti/Mommybrain/Slide_tags/EdgeR:/workspace \
#   --bind /scratch/mfafouti/Mommybrain/Slide_tags/EdgeR/out:/workspace/out \
#   --bind /scratch/mfafouti/Mommybrain/Slide_tags/EdgeR/out/new_march_26/pseudobulk_outputs:/workspace/pseudobulk_outputs \
#   --bind /scratch/mfafouti/miniforge3/envs/seurat_env:/opt/seurat_env \
#   /scratch/mfafouti/docker/edger.sif

# ========= LIBRARIES =========
library(edgeR)
library(stringr, lib.loc = "/opt/seurat_env/lib/R/library")
library(vctrs, lib.loc = "/opt/seurat_env/lib/R/library")

# ========= CONTAINER PATHS =========
input_dir <- "/workspace/out/new_march_26/pseudobulk_outputs"
output_dir <- "/workspace/out/edger_lrt"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- list(
  list(
    GroupA_name = "OIL",
    GroupA_samples =c("BC15","BC14","BC13"),
    GroupB_name = "CORT",
    GroupB_samples = c("BC28","BC3","BC9")
  )
)

# List all counts files
counts_files <- list.files(input_dir, pattern = "_counts.tsv$", full.names = TRUE)
print('The counts files found are:')
print(counts_files)

for (counts_file in counts_files) {
  celltype <- str_replace(basename(counts_file), "_counts.tsv", "")
  message("Processing cell type: ", celltype)
  
  # Load data
  counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1)
  print(colnames(counts))

   for (comp in comparisons) {
    groupA_name <- comp$GroupA_name
    groupB_name <- comp$GroupB_name
    groupA_samples <- comp$GroupA_samples
    groupB_samples <- comp$GroupB_samples

    # Keep only samples present in counts
    groupA_present <- groupA_samples[groupA_samples %in% colnames(counts)]
    groupB_present <- groupB_samples[groupB_samples %in% colnames(counts)]

    # Skip comparison if either group is empty
    if(length(groupA_present) < 2 | length(groupB_present) < 2) {
      cat("  Skipping comparison", groupA_name, "vs", groupB_name, 
          "for cell type", celltype, ": not enough samples in one of the groups.\n")
      next
    }

    # Subset counts and create group factor
    samples_sub <- c(groupA_present, groupB_present)
    counts_sub <- counts[, samples_sub, drop=FALSE]
    group <- factor(c(rep(groupA_name, length(groupA_present)),
                      rep(groupB_name, length(groupB_present))), 
                      levels = c("OIL", "CORT"))

    # --- edgeR analysis ---
    dge <- DGEList(counts = counts_sub, group = group)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~group)

    # ---- GENE LEVEL FILTERING ----
    keep <- filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]

    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, coef = 2)  # Assumes treatmentB is second level i.e CORT
  
    # Save results
    result_file <- file.path(output_dir, paste0(celltype, "_edgeR_results.tsv"))
    write.table(topTags(lrt, n = Inf)$table, result_file, sep = "\t", quote = FALSE)
   }
}

message("✅ edgeR analysis complete for all cell types.")