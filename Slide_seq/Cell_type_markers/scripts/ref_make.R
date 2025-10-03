# ========== Load libraries ==========
library(Matrix)
library(Seurat)

# ========== Set paths ==========
base_dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/Cell_type_markers/data/anndata_export"
split_dirs <- c("split1", "split2")

seurat_list <- list()

# ========== Load each split and create Seurat object ==========
for (split in split_dirs) {
  split_path <- file.path(base_dir, split)
  
  # Load expression matrix
  expr_matrix <- readMM(file.path(split_path, "matrix.mtx"))
  expr_matrix <- t(expr_matrix)  # transpose so genes are rows
  
  # Load features (genes) and barcodes (cells)
  features <- read.csv(file.path(split_path, "features.csv"), row.names = 1)
  barcodes <- read.csv(file.path(split_path, "barcodes_metadata.csv"), row.names = 1)
  
  # Set row and column names
  rownames(expr_matrix) <- rownames(features)
  colnames(expr_matrix) <- rownames(barcodes)
  
  # Create Seurat object (barcodes as metadata)
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, meta.data = barcodes)
  
  # Store in list
  seurat_list[[split]] <- seurat_obj
}

# ========== Merge the two Seurat objects ==========
merged_seurat <- merge(seurat_list[[1]], y = seurat_list[[2]], add.cell.ids = split_dirs)

# Print merged object summary
print(merged_seurat)

# ========== Save merged Seurat object ==========
saveRDS(merged_seurat, file = file.path(base_dir, "slide_tags_reference_merged.rds"))


# # PREVIEW Ref file:

# sobj <- readRDS("/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Reference_data/new/ref_slide_tags.rds")
# print(sobj)