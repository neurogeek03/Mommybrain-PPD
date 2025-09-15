# ========== Load libraries ==========
library(Matrix)
library(Seurat)
library(stringr)

# ========== Set base directory ==========
base_dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/UPDATED_anndata_extracts"
out_dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/UPDATED_mouse_query_rds"

# ========== Get all subdirectories ==========
sample_dirs <- list.dirs(base_dir, recursive = FALSE)

# ========== Loop through each sample ==========
for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  message("📦 Processing: ", sample_name)
  
  # Set file paths
  matrix_path <- file.path(sample_dir, paste0(sample_name, "_matrix.mtx"))
  genes_path <- file.path(sample_dir, paste0(sample_name, "_genes.tsv"))
  barcodes_path <- file.path(sample_dir, paste0(sample_name, "_barcodes.tsv"))
  
  # Check if all required files exist
  if (!all(file.exists(matrix_path, genes_path, barcodes_path))) {
    warning("⚠️ Missing files for ", sample_name, ", skipping...")
    next
  }

  # ========== Load expression matrix ==========
  expr_matrix <- readMM(matrix_path)
  expr_matrix <- t(expr_matrix)  # Seurat requires genes as rows

  # ========== Load gene and cell IDs ==========
  genes <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

  # ========== Assign row and column names ==========
  rownames(expr_matrix) <- genes$V1
  colnames(expr_matrix) <- barcodes$V1

  # ========== Create Seurat object ==========
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # ========== Save RDS ==========
  saveRDS(seurat_obj, file = file.path(out_dir, paste0(sample_name, ".rds")))
  message("✅ Saved: ", sample_name, ".rds")
}