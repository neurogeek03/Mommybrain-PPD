#!/usr/bin/env python3
"""
merge_and_count_subclasses.py

Load a combined Slide-seq AnnData object, apply group-specific QC filters,
and save the filtered object.

Groups & filters:
  - Neuronal     (subclass ends in 'IMN'): RCTD_singlet_score_rat > 330
  - Non-neuronal (subclass ends in 'NN') : RCTD_spot_class_rat == 'singlet'
                                           AND RCTD_singlet_score_rat > 330

Outputs (written to ../output/):
  - filtered_combined_<stamp>.h5ad     : filtered AnnData object
  - subclass_cell_counts_<stamp>.csv   : cell counts per subclass + broad class
  - umap_filtered_<stamp>.png          : UMAP coloured by subclass
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from pathlib import Path
from datetime import datetime

# ── Paths ─────────────────────────────────────────────────────────────────────
project_dir = Path(__file__).resolve().parents[1]
output_dir  = project_dir / "output"
output_dir.mkdir(exist_ok=True)

stamp       = datetime.now().strftime("%Y%m%d_%H%M%S")
output_base = output_dir

FILE_A = "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/data/All_RCTD_types_singlet_score_0_slide_seq_15.h5ad"

group_col = "RCTD_first_type_rat"

# ── Load ──────────────────────────────────────────────────────────────────────
print("Loading dataset …")
adata = sc.read_h5ad(FILE_A)
print(f"  Loaded : {adata.n_obs:>6,} cells × {adata.n_vars:>5,} genes")

# ── Split ─────────────────────────────────────────────────────────────────────
neuronal_mask = adata.obs[group_col].str.endswith(("IMN", "Glut", "Gaba"))
nn_mask       = adata.obs[group_col].str.endswith("NN")

adata_neuronal = adata[neuronal_mask].copy()
adata_nn       = adata[nn_mask].copy()

print(f"\nBefore filtering:")
print(f"  Neuronal (IMN)    : {adata_neuronal.n_obs:>6,} cells")
print(f"  Non-neuronal (NN) : {adata_nn.n_obs:>6,} cells")

# ── Filter ────────────────────────────────────────────────────────────────────
# Neuronal: RCTD_singlet_score_rat > 330
neuronal_filtered = adata_neuronal[
    adata_neuronal.obs["RCTD_singlet_score_rat"] > 330
].copy()

# Non-neuronal: RCTD_spot_class_rat == 'singlet' AND RCTD_singlet_score_rat > 330
nn_filtered = adata_nn[
    (adata_nn.obs["RCTD_spot_class_rat"] == "singlet") &
    (adata_nn.obs["RCTD_singlet_score_rat"] > 330)
].copy()

print(f"\nAfter filtering:")
print(f"  Neuronal (IMN)    : {neuronal_filtered.n_obs:>6,} cells retained")
print(f"  Non-neuronal (NN) : {nn_filtered.n_obs:>6,} cells retained")

# ── Combine ───────────────────────────────────────────────────────────────────
adata_filtered = ad.concat([neuronal_filtered, nn_filtered], join="inner", merge="same")
print(f"\n  Combined filtered : {adata_filtered.n_obs:>6,} cells × {adata_filtered.n_vars:>5,} genes")

# ── Clean var['name'] — keep gene symbol (everything before the '-') ──────────
adata_filtered.var["name"] = adata_filtered.var["name"].str.split("-").str[0]

# ── Save ──────────────────────────────────────────────────────────────────────
out_h5ad = output_base / f"filtered_combined_{stamp}.h5ad"
adata_filtered.write_h5ad(out_h5ad)
print(f"\nSaved filtered object → {out_h5ad}")

# =================== SUBCLASS EXPLORATION ===================

def get_broad_class(subclass_name):
    if subclass_name.endswith(("IMN", "Glut", "Gaba")):
        return "Neuronal"
    elif subclass_name.endswith("NN"):
        return "Non-neuronal"
    return "Unknown"

cell_counts = adata_filtered.obs[group_col].value_counts()
count_df = cell_counts.reset_index()
count_df.columns = [group_col, "n_cells"]
count_df["broad_class"] = count_df[group_col].apply(get_broad_class)
count_df = count_df.sort_values("n_cells", ascending=False)

print("\n--- Cell counts per broader class ---")
print(count_df.groupby("broad_class")["n_cells"].agg(["sum", "count"]))
print(f"\n--- All subclasses ranked by cell count (n={len(count_df)}) ---")
print(count_df.to_string(index=False))

count_df.to_csv(output_base / f"subclass_cell_counts_{stamp}.csv", index=False)
print(f"Saved subclass_cell_counts_{stamp}.csv")
