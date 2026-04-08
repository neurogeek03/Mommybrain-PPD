"""
Plot Dcx expression on UMAP with a continuous gradient.
- Zero/no-expression cells: light grey
- Dcx+ cells: yellow → red gradient (higher = darker red)

Usage:
    conda run -n anndata_env python scripts/plot_dcx_umap.py
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import os

# ── paths ──────────────────────────────────────────────────────────────────────
DATA = "data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
OUT_DIR = "out"
os.makedirs(OUT_DIR, exist_ok=True)

# ── load ───────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(DATA)
print(f"  Shape: {adata.shape}")

# ── find Dcx ──────────────────────────────────────────────────────────────────
# Locate the var column whose name contains "gene_symbol"
sym_cols = [c for c in adata.var.columns if "gene_symbol" in c.lower()]
if not sym_cols:
    raise ValueError(f"No gene_symbol column found in adata.var. Available columns: {adata.var.columns.tolist()}")
sym_col = sym_cols[0]
print(f"  Using var column: '{sym_col}'")

symbols = adata.var[sym_col].astype(str)
mask_exact = symbols.str.lower() == "dcx"
dcx_hits = adata.var_names[mask_exact].tolist()

if not dcx_hits:
    # fallback: partial match
    mask_partial = symbols.str.lower().str.contains("dcx")
    dcx_hits = adata.var_names[mask_partial].tolist()

if not dcx_hits:
    raise ValueError(f"No gene matching 'Dcx' found in adata.var['{sym_col}'].")

gene = dcx_hits[0]
gene_label = symbols[gene]  # e.g. "Dcx"
print(f"  Using gene: {gene_label} (var_name: {gene}, {len(dcx_hits)} hit(s))")

# ── check UMAP ────────────────────────────────────────────────────────────────
if "X_umap" not in adata.obsm:
    raise ValueError("X_umap not found in adata.obsm. Compute UMAP first.")

# ── extract expression ────────────────────────────────────────────────────────
import scipy.sparse as sp

idx = list(adata.var_names).index(gene)
expr = adata.X[:, idx]
if sp.issparse(expr):
    expr = np.asarray(expr.todense()).flatten()
else:
    expr = np.asarray(expr).flatten()

umap = adata.obsm["X_umap"]

# ── split zero vs expressing cells ───────────────────────────────────────────
mask_pos = expr > 0
print(f"  {gene_label}+ cells: {mask_pos.sum()} / {len(mask_pos)}")

# ── plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 7))

# background: all cells in grey
ax.scatter(
    umap[:, 0], umap[:, 1],
    s=1, c="lightgrey", linewidths=0, rasterized=True, label=f"{gene_label}\u2013"
)

# foreground: Dcx+ cells with gradient
cmap = plt.cm.YlOrRd
sc_pos = ax.scatter(
    umap[mask_pos, 0], umap[mask_pos, 1],
    s=4, c=expr[mask_pos], cmap=cmap,
    vmin=expr[mask_pos].min(), vmax=expr[mask_pos].max(),
    linewidths=0, rasterized=True, label=f"{gene_label}+"
)

cbar = plt.colorbar(sc_pos, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label(f"{gene_label} expression", fontsize=11)

ax.set_xlabel("UMAP 1", fontsize=12)
ax.set_ylabel("UMAP 2", fontsize=12)
ax.set_title(f"{gene_label} expression on UMAP", fontsize=13)
ax.axis("off")

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "dcx_umap.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path}")
plt.close()
