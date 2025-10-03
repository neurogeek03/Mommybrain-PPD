# ========== IMPORTS ==========
import scanpy as sc
import os
import pandas as pd
from scipy import io, sparse
from scipy.sparse import csr_matrix

# ========== PATHS ==========
filtered_ad = "/scratch/mfafouti/Mommybrain/Slide_seq/Cell_type_markers/data/filtered_10_subclass_sample_umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad"
base_out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Cell_type_markers/data/anndata_export"

# Load full AnnData
adata = sc.read_h5ad(filtered_ad)
print("Loaded filtered .h5ad file.")

adata.obs_names_make_unique()

# Decide split: here we split cells into first half / second half
n_cells = adata.n_obs
split_idx = n_cells // 2

adata1 = adata[:split_idx, :].copy()
adata2 = adata[split_idx:, :].copy()

# ========== HELPER FUNCTION ==========
def export_anndata(ad, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # Get expression matrix
    if ad.raw is not None:
        expr = ad.raw.X
        print("Using ad.raw.X for export.")
    elif "counts" in ad.layers:
        expr = ad.layers["counts"]
        print("Using ad.layers['counts'] for export.")
    else:
        expr = ad.X
        print("Using ad.X (ensure it is sparse).")

    # Ensure sparse
    if not sparse.issparse(expr):
        expr = csr_matrix(expr)

    # Save matrix
    io.mmwrite(os.path.join(out_dir, "matrix.mtx"), expr)
    print(f"Exported expression matrix to {os.path.join(out_dir, 'matrix.mtx')}")

    # Save features/genes
    ad.var.to_csv(os.path.join(out_dir, "features.csv"))
    print(f"Exported genes to {os.path.join(out_dir, 'features.csv')}")

    # Save barcodes
    ad.obs.to_csv(os.path.join(out_dir, "barcodes_metadata.csv"))
    print(f"Exported cell barcodes to {os.path.join(out_dir, 'barcodes_metadata.csv')}")

    # Save one metadata column
    if "RCTD_first_type_mouse" in ad.obs.columns:
        ad.obs[["RCTD_first_type_mouse"]].to_csv(
            os.path.join(out_dir, "subclass_labels.csv"),
            index=True,
            index_label="cell_label"
        )
        print(f"Exported cell metadata to {os.path.join(out_dir, 'subclass_labels.csv')}")

# ========== EXPORT BOTH SPLITS ==========
export_anndata(adata1, os.path.join(base_out_dir, "split1"))
export_anndata(adata2, os.path.join(base_out_dir, "split2"))

print("Finished exporting both sub-files.")
