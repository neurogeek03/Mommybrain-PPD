library(Seurat)
library(Matrix)

print("libaries loaded!")

# Set your path
dir <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/WMB_ref/min20"

# Define base name (must match the filenames)
basename <- "WMB"

# Load expression matrix
counts <- readMM(file.path(dir, paste0(basename, "_matrix.mtx")))

# FIX: transposing matrix
counts <- t(counts)

print("exp matrix loaded!")

# Load gene and cell names
genes <- read.delim(file.path(dir, paste0(basename, "_genes.tsv")), header = FALSE)[,1]
barcodes <- read.delim(file.path(dir, paste0(basename, "_barcodes.tsv")), header = FALSE)[,1]

# Add row and column names to matrix
rownames(counts) <- genes
colnames(counts) <- barcodes

print("row & colnames added !")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = basename)

print("sobj created!")

# Optional: Add metadata (e.g., class)
meta_file <- "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/WMB_ref/min20/WMB_class.csv"
if (file.exists(meta_file)) {
  class_meta <- read.csv(meta_file, row.names = 1)
  
  # Ensure the order of barcodes matches
  class_meta <- class_meta[colnames(seurat_obj), , drop = FALSE]
  
  # Add metadata to Seurat object
  seurat_obj <- AddMetaData(seurat_obj, metadata = class_meta)
}
  print("✅ class metadata added to Seurat object.")
} else {
  print("⚠️ class metadata file not found. Skipping metadata.")
}

# Save as RDS
saveRDS(seurat_obj, file = file.path(dir, paste0(basename,"class_meta", ".rds")))
print(paste("✅ Saved Seurat object to", file.path(dir, paste0(basename,"class_meta", ".rds"))))
