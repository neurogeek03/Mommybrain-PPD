import scanpy as sc
import pandas as pd
from pathlib import Path

# =================== INPUT ===================
h5ad_path = Path('/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/out/semi_filtered_neurons/181663_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad')

adata = sc.read_h5ad(h5ad_path)
print(f"Loaded: {h5ad_path.name}")
print(f"Total cells: {adata.n_obs}")
print(f"\nObs columns: {list(adata.obs.columns)}\n")

# =================== IDENTIFY SAMPLE COLUMN ===================
# Common sample ID column names — adjust if needed
sample_col_candidates = ['sample', 'sample_id', 'batch', 'donor', 'subject', 'orig.ident', 'library_id']

sample_col = None
for col in sample_col_candidates:
    if col in adata.obs.columns:
        sample_col = col
        break

if sample_col is None:
    print("Could not auto-detect sample column. Showing value_counts for all object/category columns:\n")
    for col in adata.obs.columns:
        if adata.obs[col].dtype in ['object', 'category']:
            n_unique = adata.obs[col].nunique()
            if 1 < n_unique <= 100:
                print(f"--- {col} ({n_unique} unique values) ---")
                print(adata.obs[col].value_counts().to_string())
                print()
else:
    print(f"Sample column detected: '{sample_col}'\n")
    counts = adata.obs[sample_col].value_counts().sort_index()
    counts.name = 'n_cells'
    counts.index.name = sample_col
    print(counts.to_string())
    print(f"\nTotal: {counts.sum()} cells across {len(counts)} samples")

    # Save to CSV
    out_csv = h5ad_path.parent / f"{h5ad_path.stem}_cell_counts_per_sample.csv"
    counts.reset_index().to_csv(out_csv, index=False)
    print(f"\nSaved to: {out_csv}")
