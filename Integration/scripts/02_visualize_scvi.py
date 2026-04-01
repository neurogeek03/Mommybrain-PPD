"""
02_visualize_scvi.py
--------------------
Visualize the scVI latent space from a saved model.

Loads merged_raw.h5ad + saved scVI model, extracts the latent
representation, computes UMAP, and generates:
  - Static PNGs (method, cell type, sample, treatment, day)
  - Interactive HTML (method, cell type) — subsampled for performance

Usage:
  python 02_visualize_scvi.py [config_path]
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.express as px
import yaml
from pathlib import Path

# =============================================================================
# CONFIG
# =============================================================================

import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    "config",
    nargs="?",
    default=Path(__file__).parent.parent / "config" / "01_scanvi.conf",
)
args = parser.parse_args()

with open(args.config) as f:
    cfg = yaml.safe_load(f)

MERGED_PATH  = cfg["merged_path"]
OUT_DIR      = Path(cfg["out_dir"])
CONFIG_NAME  = cfg["config_name"]
SCVI_PARAMS  = cfg["scvi"]
COLOR_CSV    = "/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv"

CONFIG_DIR   = OUT_DIR / CONFIG_NAME
SCVI_PATH    = CONFIG_DIR / "scvi_model"
VIZ_DIR      = CONFIG_DIR / "viz_scvi"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# LOAD DATA + MODEL
# =============================================================================

print(f"[{CONFIG_NAME}] Loading merged object...")
adata = sc.read_h5ad(MERGED_PATH)
print(f"Shape: {adata.shape}")

print("Loading saved scVI model...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer=SCVI_PARAMS["layer"],
    batch_key=SCVI_PARAMS["batch_key"],
)
scvi_model = scvi.model.SCVI.load(str(SCVI_PATH), adata=adata)
print("Model loaded.")

# =============================================================================
# LATENT REPRESENTATION + UMAP
# =============================================================================

print("Extracting scVI latent representation...")
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()

print("Computing neighbors + UMAP...")
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

# =============================================================================
# COLOR MAP (Allen Brain Cell Atlas)
# =============================================================================

color_df = pd.read_csv(COLOR_CSV, usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

def build_color_map(labels):
    cm = {}
    for label in labels:
        hex_val = label_to_hex.get(label, "999999")
        cm[label] = f"#{hex_val}" if not hex_val.startswith("#") else hex_val
    return cm

cell_type_colors = build_color_map(adata.obs["cell_type"].unique())

# =============================================================================
# STATIC PLOTS (PNG)
# =============================================================================

print("Generating static UMAP plots...")

# 1 — method + cell type side-by-side
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata, color="method", ax=axes[0], show=False,
           title="scVI latent — method")
sc.pl.umap(adata, color="cell_type", ax=axes[1], show=False,
           title="scVI latent — cell type",
           legend_loc="on data", legend_fontsize=4)
fig.suptitle(f"scVI integration check — {CONFIG_NAME}", fontsize=13)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scvi_umap_method_celltype.png", dpi=150)
plt.close(fig)
print("Saved: scvi_umap_method_celltype.png")

# 2 — treatment + day side-by-side
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata, color="treatment", ax=axes[0], show=False, title="Treatment")
sc.pl.umap(adata, color="day",       ax=axes[1], show=False, title="Day")
fig.tight_layout()
fig.savefig(VIZ_DIR / "scvi_umap_treatment_day.png", dpi=150)
plt.close(fig)
print("Saved: scvi_umap_treatment_day.png")

# 3 — per-sample (to check for sample-level clustering / batch effects)
fig, ax = plt.subplots(figsize=(12, 7))
sc.pl.umap(adata, color="sample_id", ax=ax, show=False,
           title="scVI latent — sample", legend_fontsize=6)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scvi_umap_sample.png", dpi=150)
plt.close(fig)
print("Saved: scvi_umap_sample.png")

# =============================================================================
# INTERACTIVE PLOTS (HTML) — subsampled
# =============================================================================

print("Generating interactive UMAP plots (subsampled)...")

N_SUBSAMPLE = 60_000
rng = np.random.default_rng(42)
idx = rng.choice(adata.n_obs, size=min(N_SUBSAMPLE, adata.n_obs), replace=False)
sub = adata[idx]

umap_coords = sub.obsm["X_umap"]
plot_df = pd.DataFrame({
    "UMAP1":      umap_coords[:, 0],
    "UMAP2":      umap_coords[:, 1],
    "method":     sub.obs["method"].values,
    "cell_type":  sub.obs["cell_type"].values,
    "sample_id":  sub.obs["sample_id"].values,
    "treatment":  sub.obs["treatment"].values,
    "day":        sub.obs["day"].values,
})

# Method
fig_method = px.scatter(
    plot_df, x="UMAP1", y="UMAP2",
    color="method",
    color_discrete_map={"slide_tags": "#E64B35", "slide_seq": "#4DBBD5"},
    hover_data=["cell_type", "sample_id", "treatment", "day"],
    title=f"scVI UMAP — method ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_method.update_traces(marker_size=2)
fig_method.write_html(str(VIZ_DIR / "scvi_umap_method.html"))
print("Saved: scvi_umap_method.html")

# Cell type — order by numeric prefix in label (e.g. "006_L4_5_IT_CTX_Glut" → 6)
def order_by_label_number(df, color_col):
    df = df.copy()
    df["_label_number"] = df[color_col].str.split("_").str[0].astype(int, errors="ignore")
    ordered = df.groupby(color_col, observed=True)["_label_number"].first().sort_values().index.tolist()
    df[color_col] = pd.Categorical(df[color_col], categories=ordered, ordered=True)
    df = df.drop(columns="_label_number")
    return df, ordered

plot_df_ct, ordered_labels = order_by_label_number(plot_df, "cell_type")
sub_colors = build_color_map(ordered_labels)
fig_ct = px.scatter(
    plot_df_ct, x="UMAP1", y="UMAP2",
    color="cell_type",
    color_discrete_map=sub_colors,
    category_orders={"cell_type": ordered_labels},
    hover_data=["method", "sample_id", "treatment", "day"],
    title=f"scVI UMAP — cell type ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_ct.update_traces(marker_size=2)
fig_ct.update_layout(legend=dict(font=dict(size=8), itemsizing="constant"))
fig_ct.write_html(str(VIZ_DIR / "scvi_umap_celltype.html"))
print("Saved: scvi_umap_celltype.html")

# =============================================================================
# WHAT TO LOOK FOR
# =============================================================================

print(f"""
[02_visualize_scvi.py — {CONFIG_NAME}] Done.
Outputs in: {VIZ_DIR}

What to check:
  method plot     — slide-tags and slide-seq should be MIXED, not separated
  cell type plot  — same cell types from both methods should cluster together
  sample plot     — no dominant per-sample islands (would indicate batch effect)
  treatment/day   — should NOT drive separation at this stage (model corrects for method, not biology)
""")
