"""
For each sample, count Dcx+ and Dcx- cells within the '038 DG' subclass.

Usage:
    conda run -n anndata_env python scripts/diagnostic/count_dcx_038dg.py \
        data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

SUBCLASS_COL = "subclass_name"
SAMPLE_COL   = "sample"
TARGET_CLASS = "038 DG"

parser = argparse.ArgumentParser()
parser.add_argument("h5ad", help="Path to the h5ad file")
args = parser.parse_args()

# ── load ───────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(args.h5ad)
print(f"  Shape: {adata.shape}")

# ── subset to 038 DG ──────────────────────────────────────────────────────────
if SUBCLASS_COL not in adata.obs.columns:
    raise ValueError(f"Column '{SUBCLASS_COL}' not in adata.obs. Available: {adata.obs.columns.tolist()}")

mask_class = adata.obs[SUBCLASS_COL].astype(str) == TARGET_CLASS
sub = adata[mask_class].copy()
print(f"  '{TARGET_CLASS}' cells: {sub.n_obs}")

if sub.n_obs == 0:
    vals = adata.obs[SUBCLASS_COL].unique().tolist()
    raise ValueError(f"No cells found for '{TARGET_CLASS}'. Available subclasses: {vals}")

# ── resolve Dcx ───────────────────────────────────────────────────────────────
sym_cols = [c for c in sub.var.columns if "gene_symbol" in c.lower()]
if not sym_cols:
    raise ValueError(f"No gene_symbol column in adata.var. Available: {sub.var.columns.tolist()}")
sym_col = sym_cols[0]

symbols = sub.var[sym_col].astype(str)
exact = symbols.str.lower() == "dcx"
hits = sub.var_names[exact].tolist()
if not hits:
    partial = symbols.str.lower().str.contains("dcx", na=False)
    hits = sub.var_names[partial].tolist()
if not hits:
    raise ValueError(f"Gene 'Dcx' not found in adata.var['{sym_col}'].")
if len(hits) > 1:
    print(f"  WARNING: 'Dcx' matched {len(hits)} entries: {hits} — using first.")

dcx_var = hits[0]
print(f"  Dcx var_name: {dcx_var}  symbol: {symbols[dcx_var]}")

# ── extract expression ────────────────────────────────────────────────────────
idx = list(sub.var_names).index(dcx_var)
expr = sub.X[:, idx]
if sp.issparse(expr):
    expr = np.asarray(expr.todense()).flatten()
else:
    expr = np.asarray(expr).flatten()

dcx_pos = expr > 0

# ── count per sample ──────────────────────────────────────────────────────────
if SAMPLE_COL not in sub.obs.columns:
    raise ValueError(f"Column '{SAMPLE_COL}' not in adata.obs. Available: {sub.obs.columns.tolist()}")

obs = sub.obs[[SAMPLE_COL]].copy()
obs["Dcx_status"] = np.where(dcx_pos, "Dcx+", "Dcx-")

counts = (
    obs.groupby([SAMPLE_COL, "Dcx_status"], observed=True)
    .size()
    .unstack(fill_value=0)
    .rename_axis(None, axis=1)
)

# ensure both columns exist
for col in ["Dcx+", "Dcx-"]:
    if col not in counts.columns:
        counts[col] = 0

counts = counts[["Dcx+", "Dcx-"]]
counts["total"] = counts["Dcx+"] + counts["Dcx-"]
counts["pct_Dcx+"] = (counts["Dcx+"] / counts["total"] * 100).round(1)

print(f"\nDcx+ / Dcx- breakdown per sample in '{TARGET_CLASS}':")
print(counts.to_string())