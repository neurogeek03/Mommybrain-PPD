"""
check_celltype_counts.py
------------------------
Check cell type count distributions per method to inform
filtering thresholds before scANVI training.

Usage:
  python check_celltype_counts.py
"""

import pandas as pd
import scanpy as sc
from pathlib import Path

MERGED_PATH = Path("/scratch/mfafouti/Mommybrain/Integration/out/merged_raw.h5ad")

print("Loading merged object (obs only)...")
adata = sc.read_h5ad(MERGED_PATH, backed="r")
obs = adata.obs[["method", "cell_type"]].copy()
adata.file.close()

for method in ["slide_tags", "slide_seq"]:
    counts = (
        obs[obs["method"] == method]["cell_type"]
        .value_counts()
        .rename("n_cells")
        .reset_index()
        .rename(columns={"index": "cell_type"})
    )

    print(f"\n{'='*55}")
    print(f"  {method}  —  {len(counts)} cell types, {counts['n_cells'].sum():,} cells total")
    print(f"{'='*55}")
    print(f"  {'Percentile':<15} {'n_cells'}")
    for p in [1, 5, 10, 25, 50, 75, 90, 95, 99, 100]:
        print(f"  {p:<15} {counts['n_cells'].quantile(p/100):.0f}")

    print(f"\n  Types below threshold:")
    print(f"  {'Threshold':<12} {'N types':<10} {'N cells lost':<15} {'% cells lost'}")
    for t in [3, 10, 20, 25, 50, 100]:
        below = counts[counts["n_cells"] < t]
        n_types = len(below)
        n_cells = below["n_cells"].sum()
        pct = 100 * n_cells / counts["n_cells"].sum()
        print(f"  {t:<12} {n_types:<10} {n_cells:<15,} {pct:.2f}%")

    print(f"\n  Bottom 20 cell types:")
    print(counts.tail(20).to_string(index=False))
