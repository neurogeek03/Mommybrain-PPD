"""
02_visualize_scvi.py
--------------------
Visualize the scANVI latent space from a saved integrated object.

Loads integrated.h5ad (output of 01_scanvi.py), which already contains
X_scANVI, X_umap, scanvi_prediction, and scanvi_max_prob, and generates:
  - Static PNGs (method, cell type, sample, treatment, day,
                 scanvi_prediction, prediction confidence)
  - Interactive HTML (method, cell type, scanvi_prediction) — subsampled

Usage:
  python 02_visualize_scvi.py [config_path]
"""

import numpy as np
import pandas as pd
import scanpy as sc
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

OUT_DIR      = Path(cfg["out_dir"])
CONFIG_NAME  = cfg["config_name"]
COLOR_CSV    = "/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv"

CONFIG_DIR      = OUT_DIR / CONFIG_NAME
INTEGRATED_PATH = CONFIG_DIR / "integrated.h5ad"
VIZ_DIR         = CONFIG_DIR / "viz_scanvi"
VIZ_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# LOAD DATA
# =============================================================================

print(f"[{CONFIG_NAME}] Loading integrated object...")
adata = sc.read_h5ad(INTEGRATED_PATH)
print(f"Shape: {adata.shape}")
print(f"  obsm keys: {list(adata.obsm.keys())}")
print(f"  obs cols: {list(adata.obs.columns)}")

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
           title="scANVI latent — method")
sc.pl.umap(adata, color="cell_type", ax=axes[1], show=False,
           title="scANVI latent — cell type",
           legend_loc="on data", legend_fontsize=4)
fig.suptitle(f"scANVI integration check — {CONFIG_NAME}", fontsize=13)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_method_celltype.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_method_celltype.png")

# 2 — treatment + day side-by-side
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata, color="treatment", ax=axes[0], show=False, title="Treatment")
sc.pl.umap(adata, color="day",       ax=axes[1], show=False, title="Day")
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_treatment_day.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_treatment_day.png")

# 3 — per-sample (to check for sample-level clustering / batch effects)
fig, ax = plt.subplots(figsize=(12, 7))
sc.pl.umap(adata, color="sample_id", ax=ax, show=False,
           title="scANVI latent — sample", legend_fontsize=6)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_sample.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_sample.png")

# 4 — scANVI prediction vs ground-truth cell_type side-by-side
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata, color="cell_type", ax=axes[0], show=False,
           title="cell_type (reference)", legend_loc="on data", legend_fontsize=4)
sc.pl.umap(adata, color="scanvi_prediction", ax=axes[1], show=False,
           title="scANVI prediction (all cells)", legend_loc="on data", legend_fontsize=4)
fig.suptitle(f"scANVI label transfer — {CONFIG_NAME}", fontsize=13)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_prediction.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_prediction.png")

# 5 — prediction confidence (scanvi_max_prob)
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata, color="scanvi_max_prob", ax=axes[0], show=False,
           title="Prediction confidence (all cells)",
           color_map="RdYlGn", vmin=0, vmax=1)
query_mask = adata.obs["method"] != cfg.get("reference_method", "slide_tags")
axes[1].hist(adata.obs.loc[query_mask, "scanvi_max_prob"], bins=50, edgecolor="black")
axes[1].set_xlabel("Max prediction probability")
axes[1].set_ylabel("Count")
axes[1].set_title("Prediction confidence — Slide-seq cells")
axes[1].axvline(0.5, color="red", linestyle="--", label="0.5 threshold")
axes[1].legend()
fig.suptitle(f"scANVI prediction confidence — {CONFIG_NAME}", fontsize=13)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_confidence.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_confidence.png")

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
    "UMAP1":             umap_coords[:, 0],
    "UMAP2":             umap_coords[:, 1],
    "method":            sub.obs["method"].values,
    "cell_type":         sub.obs["cell_type"].values,
    "scanvi_prediction": sub.obs["scanvi_prediction"].values,
    "scanvi_max_prob":   sub.obs["scanvi_max_prob"].values,
    "sample_id":         sub.obs["sample_id"].values,
    "treatment":         sub.obs["treatment"].values,
    "day":               sub.obs["day"].values,
})

# Method
fig_method = px.scatter(
    plot_df, x="UMAP1", y="UMAP2",
    color="method",
    color_discrete_map={"slide_tags": "#E64B35", "slide_seq": "#4DBBD5"},
    hover_data=["cell_type", "scanvi_prediction", "sample_id", "treatment", "day"],
    title=f"scANVI UMAP — method ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_method.update_traces(marker_size=2)
fig_method.write_html(str(VIZ_DIR / "scanvi_umap_method.html"))
print("Saved: scanvi_umap_method.html")

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
    hover_data=["method", "scanvi_prediction", "sample_id", "treatment", "day"],
    title=f"scANVI UMAP — cell type ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_ct.update_traces(marker_size=2)
fig_ct.update_layout(legend=dict(font=dict(size=8), itemsizing="constant"))
fig_ct.write_html(str(VIZ_DIR / "scanvi_umap_celltype.html"))
print("Saved: scanvi_umap_celltype.html")

# scANVI prediction
plot_df_pred, ordered_preds = order_by_label_number(plot_df, "scanvi_prediction")
pred_colors = build_color_map(ordered_preds)
fig_pred = px.scatter(
    plot_df_pred, x="UMAP1", y="UMAP2",
    color="scanvi_prediction",
    color_discrete_map=pred_colors,
    category_orders={"scanvi_prediction": ordered_preds},
    hover_data=["method", "cell_type", "scanvi_max_prob", "sample_id", "treatment", "day"],
    title=f"scANVI UMAP — predicted label ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_pred.update_traces(marker_size=2)
fig_pred.update_layout(legend=dict(font=dict(size=8), itemsizing="constant"))
fig_pred.write_html(str(VIZ_DIR / "scanvi_umap_prediction.html"))
print("Saved: scanvi_umap_prediction.html")

# =============================================================================
# FILTERED PLOTS — remove low-confidence Slide-seq cells (threshold = 0.8)
# Slide-tags cells are kept regardless (they are the labeled reference).
# =============================================================================

CONF_THRESHOLD = 0.8
ref_method = cfg.get("reference_method", "slide_tags")

keep_mask = (adata.obs["method"] == ref_method) | (adata.obs["scanvi_max_prob"] >= CONF_THRESHOLD)
adata_filt = adata[keep_mask].copy()

n_removed = adata.n_obs - adata_filt.n_obs
n_slideseq_total = (adata.obs["method"] != ref_method).sum()
n_slideseq_removed = (query_mask & ~keep_mask).sum()
print(f"\nFiltered: removed {n_removed} low-confidence Slide-seq cells "
      f"({n_slideseq_removed}/{n_slideseq_total}, "
      f"{100*n_slideseq_removed/n_slideseq_total:.1f}%) "
      f"with max_prob < {CONF_THRESHOLD}")
print(f"Remaining cells: {adata_filt.n_obs}")

print("Generating filtered static UMAP plots...")

# method + cell type
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata_filt, color="method", ax=axes[0], show=False,
           title=f"scANVI latent — method (conf≥{CONF_THRESHOLD})")
sc.pl.umap(adata_filt, color="cell_type", ax=axes[1], show=False,
           title=f"scANVI latent — cell type (conf≥{CONF_THRESHOLD})",
           legend_loc="on data", legend_fontsize=4)
fig.suptitle(f"scANVI filtered (conf≥{CONF_THRESHOLD}) — {CONFIG_NAME}", fontsize=13)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_method_celltype_filt.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_method_celltype_filt.png")

# prediction vs cell_type
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata_filt, color="cell_type", ax=axes[0], show=False,
           title="cell_type (reference)", legend_loc="on data", legend_fontsize=4)
sc.pl.umap(adata_filt, color="scanvi_prediction", ax=axes[1], show=False,
           title=f"scANVI prediction (conf≥{CONF_THRESHOLD})",
           legend_loc="on data", legend_fontsize=4)
fig.suptitle(f"scANVI label transfer filtered — {CONFIG_NAME}", fontsize=13)
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_prediction_filt.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_prediction_filt.png")

# treatment + day
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata_filt, color="treatment", ax=axes[0], show=False, title="Treatment")
sc.pl.umap(adata_filt, color="day",       ax=axes[1], show=False, title="Day")
fig.tight_layout()
fig.savefig(VIZ_DIR / "scanvi_umap_treatment_day_filt.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_umap_treatment_day_filt.png")

# Interactive — method (filtered)
rng2 = np.random.default_rng(42)
idx_filt = rng2.choice(adata_filt.n_obs, size=min(N_SUBSAMPLE, adata_filt.n_obs), replace=False)
sub_filt = adata_filt[idx_filt]
umap_filt = sub_filt.obsm["X_umap"]
plot_df_filt = pd.DataFrame({
    "UMAP1":             umap_filt[:, 0],
    "UMAP2":             umap_filt[:, 1],
    "method":            sub_filt.obs["method"].values,
    "cell_type":         sub_filt.obs["cell_type"].values,
    "scanvi_prediction": sub_filt.obs["scanvi_prediction"].values,
    "scanvi_max_prob":   sub_filt.obs["scanvi_max_prob"].values,
    "sample_id":         sub_filt.obs["sample_id"].values,
    "treatment":         sub_filt.obs["treatment"].values,
    "day":               sub_filt.obs["day"].values,
})

fig_filt_method = px.scatter(
    plot_df_filt, x="UMAP1", y="UMAP2",
    color="method",
    color_discrete_map={"slide_tags": "#E64B35", "slide_seq": "#4DBBD5"},
    hover_data=["cell_type", "scanvi_prediction", "sample_id", "treatment", "day"],
    title=f"scANVI UMAP — method, conf≥{CONF_THRESHOLD} ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_filt_method.update_traces(marker_size=2)
fig_filt_method.write_html(str(VIZ_DIR / "scanvi_umap_method_filt.html"))
print("Saved: scanvi_umap_method_filt.html")

plot_df_filt_ct, ordered_labels_filt = order_by_label_number(plot_df_filt, "cell_type")
filt_colors = build_color_map(ordered_labels_filt)
fig_filt_ct = px.scatter(
    plot_df_filt_ct, x="UMAP1", y="UMAP2",
    color="cell_type",
    color_discrete_map=filt_colors,
    category_orders={"cell_type": ordered_labels_filt},
    hover_data=["method", "scanvi_prediction", "sample_id", "treatment", "day"],
    title=f"scANVI UMAP — cell type, conf≥{CONF_THRESHOLD} ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_filt_ct.update_traces(marker_size=2)
fig_filt_ct.update_layout(legend=dict(font=dict(size=8), itemsizing="constant"))
fig_filt_ct.write_html(str(VIZ_DIR / "scanvi_umap_celltype_filt.html"))
print("Saved: scanvi_umap_celltype_filt.html")

plot_df_filt_pred, ordered_preds_filt = order_by_label_number(plot_df_filt, "scanvi_prediction")
filt_pred_colors = build_color_map(ordered_preds_filt)
fig_filt_pred = px.scatter(
    plot_df_filt_pred, x="UMAP1", y="UMAP2",
    color="scanvi_prediction",
    color_discrete_map=filt_pred_colors,
    category_orders={"scanvi_prediction": ordered_preds_filt},
    hover_data=["method", "cell_type", "scanvi_max_prob", "sample_id", "treatment", "day"],
    title=f"scANVI UMAP — predicted label, conf≥{CONF_THRESHOLD} ({CONFIG_NAME})",
    opacity=0.4, width=1100, height=750,
)
fig_filt_pred.update_traces(marker_size=2)
fig_filt_pred.update_layout(legend=dict(font=dict(size=8), itemsizing="constant"))
fig_filt_pred.write_html(str(VIZ_DIR / "scanvi_umap_prediction_filt.html"))
print("Saved: scanvi_umap_prediction_filt.html")

# =============================================================================
# WHAT TO LOOK FOR
# =============================================================================

print(f"""
[02_visualize_scvi.py — {CONFIG_NAME}] Done.
Outputs in: {VIZ_DIR}

What to check:
  method plot        — slide-tags and slide-seq should be MIXED, not separated
  cell type plot     — same cell types from both methods should cluster together
  prediction plot    — slide-seq predictions should match reference cell_type clusters
  confidence plot    — low-confidence (<0.5) slide-seq cells may need review
  *_filt plots       — same but with Slide-seq cells conf<{CONF_THRESHOLD} removed
  sample plot        — no dominant per-sample islands (would indicate batch effect)
  treatment/day      — should NOT drive separation at this stage
""")