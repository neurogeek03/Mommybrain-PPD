# ========== Load libraries ==========
library(Matrix)
library(Seurat)

# ========== Set paths ==========
out_dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/Slide_tags_reference"
matrix_path <- file.path(out_dir, "matrix.mtx")
features_path <- file.path(out_dir, "features.csv")
barcodes_path <- file.path(out_dir, "barcodes_metadata.csv")
mapmycells_path <- file.path(out_dir, "subclass_labels.csv")

# ========== Load components ==========

# Load expression matrix
expr_matrix <- readMM(matrix_path)
expr_matrix <- t(expr_matrix) #TRANSPOSE matrix, as seurat requires features as rows

# Load metadata
features <- read.csv(features_path, row.names = 1)
barcodes <- read.csv(barcodes_path, row.names = 1)
meta_columns <- read.csv(mapmycells_path, row.names = 1)

# Set row and column names on expression matrix
rownames(expr_matrix) <- rownames(features)  # gene IDs (Ensembl)
colnames(expr_matrix) <- rownames(barcodes)  # cell barcodes

# ========== Create Seurat object ==========
seurat_obj <- CreateSeuratObject(counts = expr_matrix, meta.data = barcodes)

seurat_obj <- AddMetaData(seurat_obj, metadata = meta_columns)

# ========== Save RDS ==========
saveRDS(seurat_obj, file = file.path(out_dir, "slide_tags_reference_unfiltered.rds"))


# # PREVIEW Ref file:

# sobj <- readRDS("/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Reference_data/new/ref_slide_tags.rds")
# print(sobj)