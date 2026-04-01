"""
00_prepare_data.py
------------------
Phase 0: data harmonization and pre-integration QC.

Run this script TWICE:
  Pass 1 (inspect_only: true):  loads data, prints obs columns and raw-count check,
                                 then exits. Use output to fill in column_map in config.
  Pass 2 (inspect_only: false): full harmonization, merge, HVG selection, QC outputs.

Usage:
  python 00_prepare_data.py [config_path]

Config is loaded from ../config/00_prepare_data.conf by default.
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from pathlib import Path

# =============================================================================
# CONFIG — loaded from YAML
# =============================================================================

parser = argparse.ArgumentParser()
parser.add_argument(
    "config",
    nargs="?",
    default=Path(__file__).parent.parent / "config" / "00_prepare_data.conf",
    help="Path to YAML config file (default: ../config/00_prepare_data.conf)",
)
args = parser.parse_args()

with open(args.config) as f:
    cfg = yaml.safe_load(f)

INSPECT_ONLY    = cfg["inspect_only"]
PATH_SLIDETAGS  = cfg["paths"]["slidetags"]
PATH_SLIDESEQ_N = cfg["paths"]["slideseq_n"]
PATH_SLIDESEQ_NN= cfg["paths"]["slideseq_nn"]
OUT_DIR         = Path(cfg["out_dir"])
EXCLUDE_SAMPLES = cfg["exclude_samples"]
N_HVG_TARGET    = cfg["n_hvg_target"]

# column_map: {target: [slidetags_col, slideseq_col]}  (null -> None)
COLUMN_MAP = {k: tuple(v) for k, v in cfg["column_map"].items()}

LABEL_MAP_SLIDETAGS = cfg.get("label_map_slidetags") or {}
LABEL_MAP_SLIDESEQ  = cfg.get("label_map_slideseq")  or {}
DEFAULTS            = cfg.get("defaults") or {}
DOUBLET_FILTER      = cfg.get("doublet_filter")

OUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# HELPERS
# =============================================================================

def check_raw(adata, name):
    """Report whether .X appears to contain raw integer counts."""
    X = adata.X
    sample = X[:200, :200] if not issparse(X) else X[:200, :200].toarray()
    is_integer = np.allclose(sample, np.round(sample))
    has_negatives = (sample < 0).any()
    print(f"\n[{name}] Raw count check:")
    print(f"  .X dtype        : {X.dtype}")
    print(f"  Values integer  : {is_integer}")
    print(f"  Has negatives   : {has_negatives}")
    print(f"  Max value (200x200 sample): {sample.max():.1f}")
    if is_integer and not has_negatives:
        print(f"  -> Looks like raw counts. OK.")
    else:
        print(f"  -> WARNING: .X does not look like raw integer counts.")
        print(f"     Check .layers for a 'counts' or 'raw_counts' layer.")
        layers = list(adata.layers.keys())
        print(f"     Available layers: {layers}")


def inspect(adata, name):
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    print(f"  Shape      : {adata.shape}")
    print(f"  obs cols   : {adata.obs.columns.tolist()}")
    print(f"  obsm keys  : {list(adata.obsm.keys())}")
    print(f"  layers     : {list(adata.layers.keys())}")
    check_raw(adata, name)


# =============================================================================
# PASS 1 — INSPECT ONLY
# =============================================================================

print("\nLoading data...")
st  = sc.read_h5ad(PATH_SLIDETAGS)
ssn = sc.read_h5ad(PATH_SLIDESEQ_N)
ssnn = sc.read_h5ad(PATH_SLIDESEQ_NN)

inspect(st,   "Slide-tags")
inspect(ssn,  "Slide-seq neurons")
inspect(ssnn, "Slide-seq non-neurons")

# Print unique values for likely metadata columns (first 20 unique values each)
print("\n--- Slide-tags: obs sample values ---")
for col in st.obs.columns:
    n_unique = st.obs[col].nunique()
    if n_unique <= 30:
        print(f"  {col} ({n_unique} unique): {sorted(st.obs[col].astype(str).unique().tolist())}")

print("\n--- Slide-seq neurons: obs sample values ---")
for col in ssn.obs.columns:
    n_unique = ssn.obs[col].nunique()
    if n_unique <= 30:
        print(f"  {col} ({n_unique} unique): {sorted(ssn.obs[col].astype(str).unique().tolist())}")

print("\n--- Slide-seq non-neurons: obs sample values ---")
for col in ssnn.obs.columns:
    n_unique = ssnn.obs[col].nunique()
    if n_unique <= 30:
        print(f"  {col} ({n_unique} unique): {sorted(ssnn.obs[col].astype(str).unique().tolist())}")

if INSPECT_ONLY:
    print("\n[INSPECT_ONLY=True] Exiting after inspection.")
    print("Fill in COLUMN_MAP and LABEL_MAP_* above, then set INSPECT_ONLY=False to run the full pipeline.")
    import sys; sys.exit(0)

# =============================================================================
# PASS 2 — FULL PIPELINE
# =============================================================================

# --- Validate COLUMN_MAP is filled in ---
st_defaults = (DEFAULTS.get("slidetags") or {}).keys()
ss_defaults = (DEFAULTS.get("slideseq") or {}).keys()
missing = [k for k, (st_col, ss_col) in COLUMN_MAP.items()
           if st_col is None and k not in st_defaults
           and ss_col is None and k not in ss_defaults]
if missing:
    raise ValueError(
        f"COLUMN_MAP has unfilled entries with no defaults: {missing}\n"
        "Run with INSPECT_ONLY=True first, then fill in the column names."
    )

# --- 1. Filter doublets from Slide-tags ---
if DOUBLET_FILTER:
    col = DOUBLET_FILTER["column"]
    keep = DOUBLET_FILTER["keep"]
    n_before = st.n_obs
    st = st[st.obs[col] == keep].copy()
    n_removed = n_before - st.n_obs
    print(f"\nDoublet filtering (Slide-tags): {n_before} -> {st.n_obs} cells ({n_removed} doublets removed)")

# --- 2. Assign method labels and harmonize obs ---
def harmonize_obs(adata, method, col_map, label_map, defaults):
    adata = adata.copy()
    adata.obs["method"] = method
    method_defaults = defaults.get("slidetags" if method == "slide_tags" else "slideseq") or {}
    for target, (src_col, _) in col_map.items():
        if method == "slide_tags":
            src = src_col
        else:
            src = _
        if src is not None and src in adata.obs.columns:
            adata.obs[target] = adata.obs[src]
        elif target in method_defaults:
            adata.obs[target] = method_defaults[target]
        elif target not in adata.obs.columns:
            adata.obs[target] = np.nan
    # Apply label harmonization
    # Normalize labels: replace spaces and hyphens with underscores for consistency
    adata.obs["cell_type"] = adata.obs["cell_type_orig"].str.replace(r"[\s\-]+", "_", regex=True)
    if label_map:
        adata.obs["cell_type"] = adata.obs["cell_type"].map(label_map).fillna(adata.obs["cell_type"])
    return adata

st   = harmonize_obs(st,   "slide_tags", COLUMN_MAP, LABEL_MAP_SLIDETAGS, DEFAULTS)
ssn  = harmonize_obs(ssn,  "slide_seq",  COLUMN_MAP, LABEL_MAP_SLIDESEQ,  DEFAULTS)
ssnn = harmonize_obs(ssnn, "slide_seq",  COLUMN_MAP, LABEL_MAP_SLIDESEQ,  DEFAULTS)

# position column: only meaningful for Slide-seq
if "position" in st.obs.columns:
    st.obs["position"] = np.nan

# --- 3. Exclude flagged samples from Slide-seq ---
n_before = ssn.n_obs + ssnn.n_obs
mask_n  = ~ssn.obs["sample_id"].isin(EXCLUDE_SAMPLES)
mask_nn = ~ssnn.obs["sample_id"].isin(EXCLUDE_SAMPLES)
ssn  = ssn[mask_n].copy()
ssnn = ssnn[mask_nn].copy()
n_after = ssn.n_obs + ssnn.n_obs
print(f"\nExcluded samples {EXCLUDE_SAMPLES} from Slide-seq: {n_before} -> {n_after} cells")

# --- 4. Concatenate Slide-seq neurons + non-neurons ---
ss = sc.concat([ssn, ssnn], join="outer", label="neuron_class",
               keys=["neuron", "non_neuron"], index_unique="-")
print(f"Slide-seq concatenated: {ss.shape}")

# --- 5. Merge Slide-tags + Slide-seq (inner join on genes) ---
n_genes_st = st.n_vars
n_genes_ss = ss.n_vars
shared_genes = st.var_names.intersection(ss.var_names)
print(f"\nGene overlap:")
print(f"  Slide-tags genes  : {n_genes_st}")
print(f"  Slide-seq genes   : {n_genes_ss}")
print(f"  Shared genes      : {len(shared_genes)} ({100*len(shared_genes)/min(n_genes_st,n_genes_ss):.1f}% of smaller set)")

if len(shared_genes) < 5000:
    print("  WARNING: fewer than 5,000 shared genes — check gene name format (e.g. Ensembl vs symbol)")

adata = sc.concat([st, ss], join="inner", label="source",
                  keys=["slide_tags", "slide_seq"], index_unique="-")
print(f"Merged AnnData: {adata.shape}")

# --- 6. Store raw counts in .layers["counts"] ---
# If .X is already raw (confirmed in Pass 1), copy directly.
# If .X is normalized, use the appropriate layer instead.
if issparse(adata.X):
    X_sample = adata.X[:200, :200].toarray()
else:
    X_sample = adata.X[:200, :200]

if np.allclose(X_sample, np.round(X_sample)) and (X_sample >= 0).all():
    adata.layers["counts"] = adata.X.copy()
    print("\n.layers['counts'] set from .X (raw counts confirmed).")
else:
    raise ValueError(
        ".X does not appear to contain raw counts after merging.\n"
        "Identify the correct layer (from Pass 1 output) and set adata.layers['counts'] manually."
    )

# =============================================================================
# QC METRICS — Checkpoint 1 outputs
# =============================================================================

print("\n--- Cell type distribution per method ---")
ct_dist = (adata.obs.groupby(["method", "cell_type"])
           .size()
           .reset_index(name="n_cells")
           .pivot(index="cell_type", columns="method", values="n_cells")
           .fillna(0)
           .astype(int))
ct_dist["total"] = ct_dist.sum(axis=1)
ct_dist = ct_dist.sort_values("total", ascending=False)
print(ct_dist.to_string())

# Flag cell types with <50 cells in any single method
for method in ["slide_tags", "slide_seq"]:
    if method in ct_dist.columns:
        low = ct_dist[ct_dist[method] < 50].index.tolist()
        if low:
            print(f"\n  WARNING: <50 cells in {method} for: {low}")

# Save cell type distribution table
ct_dist.to_csv(OUT_DIR / "00_cell_type_distribution.csv")

# --- Library size distributions ---
sc.pp.calculate_qc_metrics(adata, inplace=True)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, metric, label in zip(axes,
                              ["log1p_total_counts", "log1p_n_genes_by_counts"],
                              ["log(total counts)", "log(n genes detected)"]):
    for method, grp in adata.obs.groupby("method"):
        ax.violinplot(grp[metric].dropna(), positions=[list(adata.obs["method"].unique()).index(method)],
                      showmedians=True)
    ax.set_xticks(range(adata.obs["method"].nunique()))
    ax.set_xticklabels(adata.obs["method"].unique(), rotation=15)
    ax.set_ylabel(label)
    ax.set_title(label)
fig.suptitle("Library size distributions by method (pre-integration)")
fig.tight_layout()
fig.savefig(OUT_DIR / "00_library_size_by_method.png", dpi=150)
plt.close(fig)
print("\nSaved: 00_library_size_by_method.png")

# --- Pre-integration PCA (baseline batch separation) ---
adata_pca = adata.copy()
sc.pp.normalize_total(adata_pca, target_sum=1e4)
sc.pp.log1p(adata_pca)
sc.pp.highly_variable_genes(adata_pca, n_top_genes=3000, batch_key="method")
sc.pp.pca(adata_pca, n_comps=30, use_highly_variable=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, color in zip(axes, ["method", "cell_type"]):
    sc.pl.pca(adata_pca, color=color, ax=ax, show=False, title=f"PCA — {color} (pre-integration)")
fig.tight_layout()
fig.savefig(OUT_DIR / "00_pca_preintegration.png", dpi=150)
plt.close(fig)
print("Saved: 00_pca_preintegration.png")

# --- HVG selection for model training ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG_TARGET, batch_key="method", subset=False)
n_hvg = adata.var["highly_variable"].sum()
print(f"\nHVGs selected: {n_hvg} (target: {N_HVG_TARGET})")
print(f"Top 20 HVGs: {adata.var_names[adata.var['highly_variable']][:20].tolist()}")

# Restore raw counts to .X (scVI/scANVI require raw counts in .X or layers["counts"])
# We keep normalized in .X for any scanpy operations but scVI will use layers["counts"]
# (scvi-tools setup_anndata reads from layer if specified)

# --- Save merged object ---
out_path = OUT_DIR / "merged_raw.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved merged AnnData: {out_path}")
print(f"  Shape: {adata.shape}")
print(f"  .layers['counts']: present")
print(f"  .var['highly_variable']: present ({n_hvg} HVGs)")

print("\n[00_prepare_data.py] Done. Review QC outputs in:", OUT_DIR)
print("Next: inspect 00_pca_preintegration.png and 00_cell_type_distribution.csv")
print("      then proceed to 01_scanvi_joint_configA.py")
