"""Step-by-step debug for UCell scoring pipeline."""

import anndata as ad
import numpy as np
from pyucell import compute_ucell_scores, smooth_knn_scores

H5AD = "Oligo_states_old/data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
CELL_TYPE_COL = "subclass_name"

OLIGO_SIG = {
    "Oligo_CORT": [
        "Fkbp5", "Sgk1", "Klf9", "Ddit4", "Zbtb16",
        "Opalin-", "Elovl6-",
    ]
}

print("=== Step 1: Load full adata ===")
adata = ad.read_h5ad(H5AD)
print(f"shape: {adata.shape}")
print(f"adata.X type: {type(adata.X)}")
print(f"adata.X dtype: {adata.X.dtype}")

print("\n=== Step 2: Subset to Oligo NN ===")
mask = adata.obs[CELL_TYPE_COL] == "327 Oligo NN"
sub = adata[mask].copy()
print(f"sub shape: {sub.shape}")
print(f"sub.X type: {type(sub.X)}")
print(f"sub.X dtype: {sub.X.dtype}")
print(f"sub.X writable: {sub.X.flags.writeable if isinstance(sub.X, np.ndarray) else 'sparse'}")

print("\n=== Step 3: Check genes present in signature ===")
for g in ["Fkbp5", "Sgk1", "Klf9", "Opalin", "Elovl6"]:
    present = g in sub.var_names
    print(f"  {g}: {'FOUND' if present else 'MISSING'}")

print("\n=== Step 4: Score (n_jobs=1 to avoid forking issues) ===")
compute_ucell_scores(sub, OLIGO_SIG, n_jobs=1)
score_col = "Oligo_CORT_UCell"
print(f"Score column present: {score_col in sub.obs.columns}")
if score_col in sub.obs.columns:
    s = sub.obs[score_col]
    print(f"  n non-null: {s.notna().sum()} / {len(s)}")
    print(f"  min={s.min():.4f}  max={s.max():.4f}  mean={s.mean():.4f}")

print("\n=== Step 5: Smooth KNN scores ===")
smooth_knn_scores(sub, [score_col])
s = sub.obs[score_col]
print(f"After smooth — n non-null: {s.notna().sum()} / {len(s)}")
print(f"  min={s.min():.4f}  max={s.max():.4f}  mean={s.mean():.4f}")

print("\n=== Step 6: Assign back to adata.obs ===")
import pandas as pd
adata.obs[score_col] = np.nan
adata.obs.loc[sub.obs_names, score_col] = sub.obs[score_col].values
assigned = adata.obs[score_col]
print(f"n non-null in adata.obs: {assigned.notna().sum()} / {len(assigned)}")
print(f"  min={assigned.min():.4f}  max={assigned.max():.4f}")
