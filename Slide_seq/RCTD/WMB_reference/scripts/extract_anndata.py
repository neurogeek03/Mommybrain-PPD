# ========== IMPORTS ==========
import scanpy as sc
import os
import anndata as ad
import pandas as pd
from scipy import io, sparse

# ========== PATHS ==========
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/Slide_tags_reference"
os.makedirs(out_dir, exist_ok=True)  # Ensure output directory exists

project_path = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD"
filtered_ad = "/scratch/mfafouti/Mommybrain/Slide_tags/Integration/scanpy/integrated_Harmony_t.h5ad"

# ========== FILE CONVERSION ==========
adata = sc.read_h5ad(filtered_ad)
print("Loaded filtered .h5ad file.")
adata.obs_names_make_unique()
print(adata.var.columns)
print(adata.var_names[:5])  # using ensembl IDs for the FILTERED_persample_subclass_doublets_harmony.h5ad

print("Type of adata.X:", type(adata.X))
print(adata.layers.keys())  # check available layers

# Option 1: From raw
if adata.raw is not None:
    expr = adata.raw.X
    print("Using adata.raw.X for export.")
# Option 2: From layer
elif "counts" in adata.layers:
    expr = adata.layers["counts"]
    print("Using adata.layers['counts'] for export.")
else:
    expr = adata.X
    print("Using adata.X (ensure it is sparse).")

# Check and convert if not sparse
from scipy.sparse import csr_matrix
if not sparse.issparse(expr):
    expr = csr_matrix(expr)

# Save sparse matrix
io.mmwrite(os.path.join(out_dir, "matrix.mtx"), expr)

# # export expression matrix (X)
# io.mmwrite(os.path.join(out_dir, "matrix.mtx"), sparse.csr_matrix(adata.X))
# print("Exported expression matrix to matrix.mtx.")

# export genes
adata.var.to_csv(os.path.join(out_dir, "features.csv"))  # rows = genes
print("Exported genes to features.csv.")

# export cell barcodes
adata.obs.to_csv(os.path.join(out_dir, "barcodes_metadata.csv"))  # rows = cells
print("Exported cell barcodes to barcodes_metadata.csv.")

# export cell metadata
adata.obs[["MapMyCells_cell_type"]].to_csv(os.path.join(out_dir, "subclass_labels.csv"), index=True,index_label="cell_label")
print("Exported cell metadata to mapmycells.csv.")

# adata.obs[["class"]].to_csv(os.path.join(out_dir, "class_labels.csv"), index=True,index_label="cell_label")
# print("Exported cell metadata to mapmycells.csv.")