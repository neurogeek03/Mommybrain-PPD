"""
Step 1: Load Slide-tags data, filter to EdgeR cell types, compute dotplot statistics.

Output: out/marker_dotplot_data.csv
"""

import re
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
from datetime import datetime

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / 'code'))

from utils.markers import get_marker_union
from utils.data_utils import (
    get_broad_class,
    load_cell_types_from_edger,
    check_log_normalized,
    compute_dotplot_data,
)

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ 01_create_marker_csv started at {stamp} ------")

# =================== CONFIG ===================
GROUP_COL       = 'subclass_name'
GENE_SYMBOL_COL = 'gene_symbol'
EDGER_FOLDER    = Path('/scratch/mfafouti/CAN-2026-poster/data/EdgeR/OIL_vs_CORT')
DATA_PATH       = PROJECT / 'data' / 'slide_tags.h5ad'
OUT_DIR         = PROJECT / 'out'
OUT_DIR.mkdir(exist_ok=True, parents=True)

# =================== CELL TYPES FROM EDGER ===================
# Subclasses present in the dataset but absent from the EdgeR folder
EXTRA_SUBCLASSES = [
    '032 L5 NP CTX Glut',
    '014 LA-BLA-BMA-PA Glut',
]

selected = load_cell_types_from_edger(EDGER_FOLDER)
print(f"\nCell types from EdgeR ({len(selected)}):")
for s in selected:
    print(f"  {s}")

# =================== LOAD DATA ===================
print(f"\nLoading {DATA_PATH} ...")
adata = sc.read_h5ad(DATA_PATH)
print(adata)

if GENE_SYMBOL_COL in adata.var.columns:
    adata.var_names = adata.var[GENE_SYMBOL_COL].astype(str)
    adata.var_names_make_unique()
    print(f"var_names → '{GENE_SYMBOL_COL}'")
else:
    print(f"WARNING: '{GENE_SYMBOL_COL}' not in adata.var; keeping existing var_names")
    print(f"  adata.var columns: {adata.var.columns.tolist()}")

# =================== TOP-40 OVERLAP DIAGNOSTIC ===================
cell_counts = adata.obs[GROUP_COL].value_counts()
top40 = set(cell_counts.head(40).index)

def _normalise(s):
    return re.sub(r'[ /]', '_', s)

obs_labels = set(adata.obs[GROUP_COL].unique())
norm_to_obs = {_normalise(v): v for v in obs_labels}
selected_in_adata = [norm_to_obs[_normalise(s)] for s in selected if _normalise(s) in norm_to_obs]
for s in EXTRA_SUBCLASSES:
    if s in obs_labels and s not in selected_in_adata:
        selected_in_adata.append(s)

in_top40  = [s for s in selected_in_adata if s in top40]
not_top40 = [s for s in selected_in_adata if s not in top40]

print(f"\n--- Top-40 overlap ({len(in_top40)}/{len(selected_in_adata)} EdgeR subclasses in top 40) ---")
print(f"\n  IN top 40:")
for s in in_top40:
    print(f"    {cell_counts[s]:6d}  {s}")
print(f"\n  NOT in top 40:")
for s in not_top40:
    print(f"    {cell_counts[s]:6d}  {s}")

print(f"\n  Top 40 subclasses NOT in EdgeR list:")
edger_norm = {_normalise(s) for s in selected_in_adata}
for s in cell_counts.head(40).index:
    if _normalise(s) not in edger_norm:
        print(f"    {cell_counts[s]:6d}  {s}")

# =================== NORMALIZE ===================
if not check_log_normalized(adata):
    print("Normalizing and log-transforming...")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

# =================== FILTER SUBCLASSES ===================
matched = selected_in_adata
unmatched_edger = [s for s in selected if _normalise(s) not in norm_to_obs]
if unmatched_edger:
    print(f"\nWARNING: {len(unmatched_edger)} EdgeR names not matched in adata:")
    for s in unmatched_edger:
        print(f"  {s}")

adata = adata[adata.obs[GROUP_COL].isin(matched)].copy()
print(f"\nFiltered to {len(matched)} subclasses, {adata.n_obs} cells")

# =================== SUBCLASS DISPLAY ORDER ===================
broad_order = ['Glutamatergic', 'GABAergic', 'Non-neuronal', 'Unknown']
groups = adata.obs[GROUP_COL].unique().tolist()
subclass_order = sorted(groups, key=lambda s: (broad_order.index(get_broad_class(s)), s))
adata.obs[GROUP_COL] = pd.Categorical(
    adata.obs[GROUP_COL], categories=subclass_order, ordered=True
)
print(f"\nSubclass display order ({len(subclass_order)}):")
for s in subclass_order:
    print(f"  {get_broad_class(s):16s}  {s}")

# =================== MARKERS ===================
marker_genes = get_marker_union()
print(f"\nMarker union: {len(marker_genes)} genes")
print(f"  {marker_genes}")

# =================== COMPUTE & EXPORT ===================
print("\nComputing dotplot data...")
df = compute_dotplot_data(adata, marker_genes, GROUP_COL, subclass_order)

if df.empty:
    raise SystemExit("ERROR: dotplot DataFrame is empty — no subclasses matched. "
                     "Check the subclass name diagnostic above.")

out_csv = OUT_DIR / 'marker_dotplot_data.csv'
df.to_csv(out_csv, index=False)
print(f"\nSaved {out_csv}")
print(f"  {df['gene'].nunique()} genes x {df['subclass'].nunique()} subclasses  "
      f"({len(df)} rows)")
print(f"\n------ 01_create_marker_csv completed ------")
