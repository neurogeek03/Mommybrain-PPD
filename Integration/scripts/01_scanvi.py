"""
01_scanvi.py
------------
Train scVI then scANVI for joint integration of Slide-tags + Slide-seq.

Slide-tags cells provide supervised labels (cell_type);
Slide-seq cells are unlabeled ("unknown") for semi-supervised learning.

Usage:
  python 01_scanvi.py [config_path]
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# CONFIG
# =============================================================================

parser = argparse.ArgumentParser()
parser.add_argument(
    "config",
    nargs="?",
    default=Path(__file__).parent.parent / "config" / "01_scanvi.conf",
)
parser.add_argument(
    "--resume",
    action="store_true",
    help="Skip scVI training and load saved model from out_dir/config_name/scvi_model",
)
args = parser.parse_args()

with open(args.config) as f:
    cfg = yaml.safe_load(f)

MERGED_PATH      = cfg["merged_path"]
OUT_DIR          = Path(cfg["out_dir"])
CONFIG_NAME      = cfg["config_name"]
SCVI_PARAMS      = cfg["scvi"]
SCANVI_PARAMS    = cfg["scanvi"]
REF_METHOD       = cfg["reference_method"]
EVAL_METRICS     = cfg.get("eval_metrics", True)

# Config-specific output directory
CONFIG_DIR = OUT_DIR / CONFIG_NAME
CONFIG_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# LOAD DATA
# =============================================================================

print(f"\n[{CONFIG_NAME}] Loading merged object...")
adata = sc.read_h5ad(MERGED_PATH)
print(f"Shape: {adata.shape}")

# --- Set up scANVI labels ---
# Reference method (Slide-tags) keeps cell_type labels
# Query method (Slide-seq) gets "unknown" for semi-supervised learning
# Rare Slide-tags types (<min_cells_per_label) are also remapped to unlabeled
# to prevent NaN during training; cells are retained in the object.
unlabeled = SCANVI_PARAMS["unlabeled_category"]
min_cells = SCANVI_PARAMS.get("min_cells_per_label", 10)

adata.obs["scanvi_labels"] = adata.obs["cell_type"].astype(str).copy()
query_mask = adata.obs["method"] != REF_METHOD
adata.obs.loc[query_mask, "scanvi_labels"] = unlabeled

label_counts = adata.obs.loc[~query_mask, "scanvi_labels"].value_counts()
rare_labels = label_counts[label_counts < min_cells].index.tolist()
if rare_labels:
    print(f"Remapping {len(rare_labels)} rare label(s) with <{min_cells} cells to '{unlabeled}' (cells retained): {rare_labels}")
    adata.obs.loc[adata.obs["scanvi_labels"].isin(rare_labels), "scanvi_labels"] = unlabeled

n_labeled = (~query_mask).sum()
n_unlabeled = query_mask.sum()
n_labels = adata.obs.loc[~query_mask, "scanvi_labels"].nunique()
print(f"Labels: {n_labeled} labeled ({n_labels} types), {n_unlabeled} unlabeled")

# =============================================================================
# scVI
# =============================================================================

scvi_path = CONFIG_DIR / "scvi_model"

if args.resume:
    print(f"\n--- Loading saved scVI model (--resume) ---")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=SCVI_PARAMS["layer"],
        batch_key=SCVI_PARAMS["batch_key"],
    )
    scvi_model = scvi.model.SCVI.load(str(scvi_path), adata=adata)
    print(f"Loaded: {scvi_path}")
else:
    print(f"\n--- Training scVI ---")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=SCVI_PARAMS["layer"],
        batch_key=SCVI_PARAMS["batch_key"],
    )

    scvi_model = scvi.model.SCVI(
        adata,
        n_layers=SCVI_PARAMS["n_layers"],
        n_latent=SCVI_PARAMS["n_latent"],
        gene_likelihood=SCVI_PARAMS["gene_likelihood"],
    )

    print(scvi_model)
    scvi_model.train(
        max_epochs=SCVI_PARAMS["max_epochs"],
        early_stopping=SCVI_PARAMS["early_stopping"],
    )

    scvi_model.save(str(scvi_path), overwrite=True)
    print(f"scVI model saved: {scvi_path}")

    train_history = scvi_model.history
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(train_history["elbo_train"]["elbo_train"], label="train")
    if "elbo_validation" in train_history:
        ax.plot(train_history["elbo_validation"]["elbo_validation"], label="validation")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("ELBO")
    ax.set_title(f"scVI training — {CONFIG_NAME}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(CONFIG_DIR / "scvi_training_history.png", dpi=150)
    plt.close(fig)
    print("Saved: scvi_training_history.png")

# =============================================================================
# scANVI
# =============================================================================

print(f"\n--- Training scANVI ---")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    unlabeled_category=unlabeled,
    labels_key="scanvi_labels",
)

print(scanvi_model)
scanvi_train_kwargs = dict(
    max_epochs=SCANVI_PARAMS["max_epochs"],
    early_stopping=SCANVI_PARAMS["early_stopping"],
)
if "early_stopping_patience" in SCANVI_PARAMS:
    scanvi_train_kwargs["early_stopping_patience"] = SCANVI_PARAMS["early_stopping_patience"]
scanvi_model.train(**scanvi_train_kwargs)

# Save scANVI model
scanvi_path = CONFIG_DIR / "scanvi_model"
scanvi_model.save(str(scanvi_path), overwrite=True)
print(f"scANVI model saved: {scanvi_path}")

# Training history
train_history = scanvi_model.history
fig, axes = plt.subplots(1, 2, figsize=(14, 4))
for ax, key, label in zip(
    axes,
    ["elbo_train", "train_classification_loss"],
    ["ELBO", "Classification loss"],
):
    if key in train_history:
        ax.plot(train_history[key][key])
    ax.set_xlabel("Epoch")
    ax.set_ylabel(label)
    ax.set_title(f"scANVI {label} — {CONFIG_NAME}")
fig.tight_layout()
fig.savefig(CONFIG_DIR / "scanvi_training_history.png", dpi=150)
plt.close(fig)
print("Saved: scanvi_training_history.png")

# =============================================================================
# LATENT SPACE & PREDICTIONS
# =============================================================================

print(f"\n--- Extracting latent space ---")
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
adata.obs["scanvi_prediction"] = scanvi_model.predict()

# Prediction confidence for Slide-seq (unlabeled) cells
pred_probs = scanvi_model.predict(soft=True)
adata.obs["scanvi_max_prob"] = pred_probs.max(axis=1)

# Report prediction stats for query cells
query_preds = adata.obs.loc[query_mask, "scanvi_prediction"]
query_conf = adata.obs.loc[query_mask, "scanvi_max_prob"]
print(f"\nSlide-seq predictions:")
print(f"  Unique types predicted: {query_preds.nunique()}")
print(f"  Confidence: mean={query_conf.mean():.3f}, median={query_conf.median():.3f}")
print(f"  Low confidence (<0.5): {(query_conf < 0.5).sum()} ({100*(query_conf < 0.5).mean():.1f}%)")

# --- UMAP on scANVI latent ---
print("\nComputing UMAP...")
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

# --- Save integrated object ---
out_path = CONFIG_DIR / "integrated.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved: {out_path}")
print(f"  Shape: {adata.shape}")
print(f"  obsm['X_scANVI']: {adata.obsm['X_scANVI'].shape}")
print(f"  obs['scanvi_prediction']: present")

# =============================================================================
# EVALUATION (scib metrics)
# =============================================================================

if EVAL_METRICS:
    print("\n--- Computing scib metrics ---")
    try:
        import scib

        # Bio conservation
        ari = scib.metrics.ari(adata, cluster_key="scanvi_prediction", label_key="cell_type")
        nmi = scib.metrics.nmi(adata, cluster_key="scanvi_prediction", label_key="cell_type")

        # Batch correction (use scANVI embedding)
        silhouette_batch = scib.metrics.silhouette_batch(
            adata, batch_key="method", label_key="cell_type", embed="X_scANVI"
        )
        ilisi = scib.metrics.ilisi_graph(
            adata, batch_key="method", type_="embed", use_rep="X_scANVI"
        )

        metrics = {
            "config": CONFIG_NAME,
            "ARI": ari,
            "NMI": nmi,
            "silhouette_batch": silhouette_batch,
            "iLISI": ilisi,
        }
        metrics_df = pd.DataFrame([metrics])
        metrics_path = CONFIG_DIR / "scib_metrics.csv"
        metrics_df.to_csv(metrics_path, index=False)
        print(f"\nscib metrics saved: {metrics_path}")
        print(metrics_df.to_string(index=False))

    except ImportError:
        print("WARNING: scib not installed, skipping metrics. Install with: pip install scib")
    except Exception as e:
        print(f"WARNING: scib metrics failed: {e}")

# =============================================================================
# QC PLOTS
# =============================================================================

print("\n--- Generating QC plots ---")

# UMAP by method
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
sc.pl.umap(adata, color="method", ax=axes[0], show=False, title="Method")
sc.pl.umap(adata, color="scanvi_prediction", ax=axes[1], show=False,
           title="scANVI prediction", legend_loc="on data", legend_fontsize=4)
fig.suptitle(f"scANVI integration — {CONFIG_NAME}")
fig.tight_layout()
fig.savefig(CONFIG_DIR / "umap_integration.png", dpi=150)
plt.close(fig)
print("Saved: umap_integration.png")

# UMAP by treatment and day
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
sc.pl.umap(adata, color="treatment", ax=axes[0], show=False, title="Treatment")
sc.pl.umap(adata, color="day", ax=axes[1], show=False, title="Day")
fig.tight_layout()
fig.savefig(CONFIG_DIR / "umap_treatment_day.png", dpi=150)
plt.close(fig)
print("Saved: umap_treatment_day.png")

# Prediction confidence distribution
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(adata.obs.loc[query_mask, "scanvi_max_prob"], bins=50, edgecolor="black")
ax.set_xlabel("Max prediction probability")
ax.set_ylabel("Count")
ax.set_title(f"scANVI prediction confidence (Slide-seq) — {CONFIG_NAME}")
ax.axvline(0.5, color="red", linestyle="--", label="0.5 threshold")
ax.legend()
fig.tight_layout()
fig.savefig(CONFIG_DIR / "prediction_confidence.png", dpi=150)
plt.close(fig)
print("Saved: prediction_confidence.png")

print(f"\n[01_scanvi.py — {CONFIG_NAME}] Done.")
print(f"Outputs in: {CONFIG_DIR}")
