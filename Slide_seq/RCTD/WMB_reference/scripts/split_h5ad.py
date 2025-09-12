import scanpy as sc
import numpy as np
import os

# === Parameters ===
adata = sc.read_h5ad("/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference_coronal/out/1M_celltypes_WMB_high_qual_cells_all_regions_max200k.h5ad")
n_splits = 10
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference_coronal/out"
os.makedirs(out_dir, exist_ok=True)

# === Shuffle indices and split ===
np.random.seed(42)  # for reproducibility
indices = np.random.permutation(adata.n_obs)  # shuffle all cell indices
splits = np.array_split(indices, n_splits)    # split into 10 parts

# === Save each subset ===
for i, idx in enumerate(splits, start=1):
    subset = adata[idx, :].copy()
    subset.write_h5ad(os.path.join(out_dir, f"subset_{i}.h5ad"))

print(f"✅ Saved {n_splits} random subsets into {out_dir}/")