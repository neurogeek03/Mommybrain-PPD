# ── Imports ──────────────────────────────────────────────────────────────────
import scanpy as sc
import pandas as pd
import numpy as np
import scanorama
from pathlib import Path
from sklearn.metrics.pairwise import cosine_distances

# ── Paths (fill in) ──────────────────────────────────────────────────────────
REFERENCE_PATH = "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/classifier/data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"   # path to reference .h5ad
QUERY_PATH     = "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/out/ultra_filtered_NN/78689_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad"   # path to query .h5ad (non-neuronal cells)
OUTPUT_PATH    = '/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/classifier/out/78689_label_transfer_scanorama.h5ad'  # path to save output query .h5ad

# ── Column names (fill in) ───────────────────────────────────────────────────
LABEL_COL      = "subclass_name"   # obs column in reference with cell-type labels (e.g. "cell_type")
N_TOP_GENES    = 2000

# ── Load ──────────────────────────────────────────────────────────────────────
adata_ref   = sc.read_h5ad(REFERENCE_PATH)
adata_query = sc.read_h5ad(QUERY_PATH)

# ── Preprocess query only (reference is already processed) ───────────────────
adata_ref.var_names_make_unique()
sc.pp.highly_variable_genes(adata_ref, flavor="seurat", n_top_genes=N_TOP_GENES, inplace=True)

sc.pp.highly_variable_genes(adata_query, flavor="seurat", n_top_genes=N_TOP_GENES, inplace=True)

# ── Scanorama integration ─────────────────────────────────────────────────────
adatas_integrated = scanorama.correct_scanpy([adata_ref, adata_query], return_dimred=True)

# ── Cosine distances in shared embedding ──────────────────────────────────────
X_ref   = adatas_integrated[0].obsm["X_scanorama"]
X_query = adatas_integrated[1].obsm["X_scanorama"]

distances = 1 - cosine_distances(X_ref, X_query)
# shape: (n_ref_cells, n_query_cells)

# ── Label transfer ────────────────────────────────────────────────────────────
def label_transfer(dist, labels):
    """Weighted label propagation via cosine similarity."""
    lab = pd.get_dummies(labels).to_numpy().T
    class_prob = lab @ dist                            # (n_classes, n_query)
    norm = np.linalg.norm(class_prob, 2, axis=0)
    class_prob = class_prob / norm
    class_prob = (class_prob.T - class_prob.min(1)) / np.ptp(class_prob, axis=1)
    return class_prob                                  # (n_query, n_classes)

class_prob = label_transfer(distances, adata_ref.obs[LABEL_COL])

# ── Store in query obs ────────────────────────────────────────────────────────
label_categories = sorted(adata_ref.obs[LABEL_COL].cat.categories)
cp_df = pd.DataFrame(class_prob, columns=label_categories, index=adata_query.obs.index)

adata_query_out = adata_query.copy()
adata_query_out.obsm["X_scanorama"] = adatas_integrated[1].obsm["X_scanorama"]
adata_query_out.obs = pd.concat([adata_query_out.obs, cp_df], axis=1)
adata_query_out.obs["predicted_label"] = cp_df.idxmax(axis=1)

# Sanitize column names with forward slashes (not allowed in HDF5 keys)
safe_map = {col: col.replace("/", "_") for col in cp_df.columns if "/" in col}
if safe_map:
    adata_query_out.obs.rename(columns=safe_map, inplace=True)
    adata_query_out.obs["predicted_label"] = adata_query_out.obs["predicted_label"].map(
        lambda x: safe_map.get(x, x)
    )

# ── Save ──────────────────────────────────────────────────────────────────────
if '_index' in adata_query_out.obs.columns:
    adata_query_out.obs = adata_query_out.obs.drop(columns=['_index'])
adata_query_out.write_h5ad(OUTPUT_PATH)
print(f"Saved to {OUTPUT_PATH}")

# ── Visualize shared embedding ────────────────────────────────────────────────
adatas_integrated[0].obs["source"] = "reference"
adatas_integrated[1].obs["source"] = "query"

adata_combined = sc.concat([adatas_integrated[0], adatas_integrated[1]])

sc.pp.neighbors(adata_combined, use_rep="X_scanorama")
sc.tl.umap(adata_combined)

# Colored by source (reference vs query) — checks integration quality
sc.pl.umap(adata_combined, color="source", save="_scanorama_source.png")

# Reference cells colored by label; query cells grey (NaN) — checks label transfer
sc.pl.umap(adata_combined, color=LABEL_COL, save="_scanorama_labels.png")
