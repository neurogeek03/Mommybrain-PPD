"""
plot_pca_interactive.py
-----------------------
Interactive PCA plot of the merged pre-integration object using plotly.
Colors from the Allen Brain Cell Atlas annotation term CSV.

Usage:
  python plot_pca_interactive.py
"""

import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px
from pathlib import Path
from scipy.sparse import issparse

# --- Config ---
MERGED_PATH = Path("/scratch/mfafouti/Mommybrain/Integration/out/merged_raw.h5ad")
COLOR_CSV   = "/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv"
OUT_DIR     = Path("/scratch/mfafouti/Mommybrain/Integration/out")

# --- Load merged object ---
print("Loading merged object...")
adata = sc.read_h5ad(MERGED_PATH)
print(f"Shape: {adata.shape}")

# --- PCA on normalized data ---
print("Running PCA...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30, use_highly_variable=True)

# --- Build color map ---
color_df = pd.read_csv(COLOR_CSV, usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)

def get_num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except:
        return -1

color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)
color_df = color_df.sort_values("num_prefix")
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

# --- Build plot dataframe ---
pca_coords = adata.obsm["X_pca"][:, :2]
plot_df = pd.DataFrame({
    "PC1": pca_coords[:, 0],
    "PC2": pca_coords[:, 1],
    "cell_type": adata.obs["cell_type"].values,
    "method": adata.obs["method"].values,
    "sample_id": adata.obs["sample_id"].values,
    "treatment": adata.obs["treatment"].values,
    "day": adata.obs["day"].values,
})

# Subsample for performance (410k points is too heavy for plotly)
n_sample = min(50000, len(plot_df))
print(f"Subsampling to {n_sample} cells for interactive plot...")
rng = np.random.default_rng(42)
idx = rng.choice(len(plot_df), size=n_sample, replace=False)
plot_df = plot_df.iloc[idx].reset_index(drop=True)

# Map colors — fall back to grey for unmapped labels
unique_types = plot_df["cell_type"].unique()
color_map = {}
for ct in unique_types:
    if ct in label_to_hex:
        color_map[ct] = f"#{label_to_hex[ct]}" if not label_to_hex[ct].startswith("#") else label_to_hex[ct]
    else:
        color_map[ct] = "#999999"

# --- Plot: colored by cell type ---
print("Generating cell type plot...")
fig_ct = px.scatter(
    plot_df,
    x="PC1", y="PC2",
    color="cell_type",
    color_discrete_map=color_map,
    hover_data=["method", "sample_id", "treatment", "day"],
    title="PCA — cell type (pre-integration)",
    opacity=0.5,
    width=1200, height=800,
)
fig_ct.update_traces(marker_size=2)
fig_ct.update_layout(legend=dict(font=dict(size=8), itemsizing="constant"))
ct_path = OUT_DIR / "00_pca_preintegration_celltype.html"
fig_ct.write_html(str(ct_path))
print(f"Saved: {ct_path}")

# --- Plot: colored by method ---
print("Generating method plot...")
fig_method = px.scatter(
    plot_df,
    x="PC1", y="PC2",
    color="method",
    color_discrete_map={"slide_tags": "#E64B35", "slide_seq": "#4DBBD5"},
    hover_data=["cell_type", "sample_id", "treatment", "day"],
    title="PCA — method (pre-integration)",
    opacity=0.5,
    width=1200, height=800,
)
fig_method.update_traces(marker_size=2)
method_path = OUT_DIR / "00_pca_preintegration_method.html"
fig_method.write_html(str(method_path))
print(f"Saved: {method_path}")

print("Done.")
