"""
Count cells that are Dcx+ / Prox1+ / Calb1-, grouped by subclass_name.

Usage:
    conda run -n anndata_env python scripts/count_dcx_prox1_calb1.py
"""

import numpy as np
import scanpy as sc
import scipy.sparse as sp

# ── paths ──────────────────────────────────────────────────────────────────────
DATA = "data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"

QUERY_GENES = {
    "Dcx":   "positive",
    "Prox1": "positive",
    "Calb1": "negative",
}

SUBCLASS_COL = "subclass_name"

# ── load ───────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(DATA)
print(f"  Shape: {adata.shape}")

# ── locate gene_symbol column ──────────────────────────────────────────────────
sym_cols = [c for c in adata.var.columns if "gene_symbol" in c.lower()]
if not sym_cols:
    raise ValueError(f"No gene_symbol column found in adata.var. Available: {adata.var.columns.tolist()}")
sym_col = sym_cols[0]
print(f"  Using var column: '{sym_col}'")

symbols = adata.var[sym_col].astype(str)

# ── resolve each gene ──────────────────────────────────────────────────────────
def resolve_gene(name):
    exact = symbols.str.lower() == name.lower()
    hits = adata.var_names[exact].tolist()
    if not hits:
        partial = symbols.str.lower().str.contains(name.lower())
        hits = adata.var_names[partial].tolist()
    if not hits:
        raise ValueError(f"No gene matching '{name}' found in adata.var['{sym_col}'].")
    if len(hits) > 1:
        print(f"  WARNING: '{name}' matched {len(hits)} entries: {hits} — using first.")
    label = symbols[hits[0]]
    print(f"  Resolved '{name}' → var_name: {hits[0]}  symbol: {label}")
    return hits[0], label

# ── extract expression ─────────────────────────────────────────────────────────
def get_expr(var_name):
    idx = list(adata.var_names).index(var_name)
    col = adata.X[:, idx]
    return np.asarray(col.todense()).flatten() if sp.issparse(col) else np.asarray(col).flatten()

# ── case-insensitive symbol search (diagnostic) ────────────────────────────────
print("\nSearching for all symbols containing each query gene name (case-insensitive):")
for gene_name in QUERY_GENES:
    hits = adata.var[symbols.str.lower().str.contains(gene_name.lower(), na=False)][sym_col].tolist()
    print(f"  '{gene_name}' matches: {hits if hits else 'NONE'}")

# ── build filter mask ──────────────────────────────────────────────────────────
print("\nGene resolution & expression summary:")
mask = np.ones(adata.n_obs, dtype=bool)

for gene_name, direction in QUERY_GENES.items():
    var_name, label = resolve_gene(gene_name)
    expr = get_expr(var_name)
    n_expressing = (expr > 0).sum()
    print(f"    {label}: {n_expressing} / {adata.n_obs} cells express this gene")
    if direction == "positive":
        mask &= expr > 0
    else:
        mask &= expr == 0

# ── report total ───────────────────────────────────────────────────────────────
print(f"\nDcx+ / Prox1+ / Calb1-: {mask.sum()} cells  ({mask.sum() / adata.n_obs * 100:.2f}% of {adata.n_obs} total)")

# ── group by subclass_name ─────────────────────────────────────────────────────
if SUBCLASS_COL not in adata.obs.columns:
    print(f"\nWARNING: '{SUBCLASS_COL}' not found in adata.obs. Available columns: {adata.obs.columns.tolist()}")
else:
    print(f"\nBreakdown by {SUBCLASS_COL}:")
    subclass = adata.obs[SUBCLASS_COL]
    counts = subclass[mask].value_counts().sort_values(ascending=False)
    for sc_name, n in counts.items():
        pct = n / mask.sum() * 100
        print(f"  {sc_name:<40} {n:>6}  ({pct:.1f}%)")
