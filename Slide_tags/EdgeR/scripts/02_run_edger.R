# ==========================================================================
# Title:        Run DE using EdgeR
# Description:  Running EdgeR LRT model and saving the results as csv files.
#               Pipeline: DGEList → lib.size filter → filterByExpr → TMM norm
#                         → estimateDisp → glmFit → glmLRT
# Author:       Maria Eleni Fafouti
# Date:         2025-06-18
# ==========================================================================

# bash
# cd /scratch/mfafouti/Mommybrain-PPD/Slide_tags/EdgeR

# apptainer shell \
#   --bind /scratch/mfafouti/Mommybrain-PPD/Slide_tags/EdgeR:/workspace \
#   --bind /scratch/mfafouti/Mommybrain-PPD/Slide_tags/EdgeR/out:/workspace/out \
#   --bind /scratch/mfafouti/miniforge3/envs/seurat_env:/opt/seurat_env \
#   /scratch/mfafouti/docker/edger.sif

# ========= LIBRARIES =========
library(edgeR)
library(stringr, lib.loc = "/opt/seurat_env/lib/R/library")
library(vctrs, lib.loc = "/opt/seurat_env/lib/R/library")

# ========= CONTAINER PATHS =========
edger_model      <- tolower(Sys.getenv("EDGER_MODEL", unset = "lrt"))
if (!edger_model %in% c("lrt", "qlf")) stop("EDGER_MODEL must be 'lrt' or 'qlf'")
message("EdgeR model: ", edger_model)

input_dir        <- "/workspace/out/pseudobulk_outputs"
output_dir       <- "/workspace/out/edger_out"
comparisons_file <- "/workspace/comparisons.csv"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ========= LOAD COMPARISONS =========
# CSV must have columns: GroupA_name, GroupA_samples, GroupB_name, GroupB_samples
# Samples within a group are comma-separated (e.g. "BC15,BC14,BC13")
comp_table <- read.csv(comparisons_file, stringsAsFactors = FALSE)

comparisons <- lapply(seq_len(nrow(comp_table)), function(i) {
  row <- comp_table[i, ]
  list(
    GroupA_name    = trimws(row$GroupA_name),
    GroupA_samples = trimws(strsplit(row$GroupA_samples, ",")[[1]]),
    GroupB_name    = trimws(row$GroupB_name),
    GroupB_samples = trimws(strsplit(row$GroupB_samples, ",")[[1]])
  )
})

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
                      levels = c(groupA_name, groupB_name))

    tryCatch({
      # --- edgeR analysis ---
      dge <- DGEList(counts = counts_sub, group = group)

      # ---- SAMPLE FILTER: drop pseudobulk samples with very low total counts ----
      keep.samples <- dge$samples$lib.size > 5e4
      if (sum(keep.samples) < ncol(dge)) {
        dropped <- colnames(dge)[!keep.samples]
        cat("  Dropping low-lib-size samples:", paste(dropped, collapse = ", "), "\n")
        dge <- dge[, keep.samples, keep.lib.sizes = FALSE]
        group <- dge$samples$group
        groupA_present <- groupA_present[groupA_present %in% colnames(dge)]
        groupB_present <- groupB_present[groupB_present %in% colnames(dge)]
        if (length(groupA_present) < 2 | length(groupB_present) < 2) {
          cat("  Skipping after lib.size filter: insufficient samples.\n")
          return(NULL)
        }
      }

      # ---- GENE FILTER: remove lowly expressed genes on pseudobulk counts ----
      keep.genes <- filterByExpr(dge, group = group)
      cat("  Keeping", sum(keep.genes), "of", nrow(dge), "genes after filterByExpr.\n")
      dge <- dge[keep.genes, , keep.lib.sizes = FALSE]

      # ---- TMM NORMALIZATION ----
      dge <- calcNormFactors(dge)

      design <- model.matrix(~group)
      dge    <- estimateDisp(dge, design = design, robust = TRUE)

      if (edger_model == "qlf") {
        fit    <- glmQLFit(dge, design = design, robust = TRUE)
        result <- glmQLFTest(fit, coef = 2)
      } else {
        fit    <- glmFit(dge, design = design)
        result <- glmLRT(fit, coef = 2)
      }

      # Save results — one subdirectory per comparison
      comp_dir <- file.path(output_dir, paste0(groupA_name, "_vs_", groupB_name))
      dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
      result_file <- file.path(comp_dir, paste0(celltype, "_edgeR_results.tsv"))
      write.table(topTags(result, n = Inf)$table, result_file, sep = "\t", quote = FALSE)
    }, error = function(e) {
      cat("  SKIPPING", groupA_name, "vs", groupB_name, "for cell type", celltype,
          "- error:", conditionMessage(e), "\n")
    })
   }
}

message(" edgeR analysis complete for all cell types.")