import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata 
from pathlib import Path

# params
min_cells = 10
samples = 2
cutoff_bootstrap = 0.5

min_cells_DE = 10

# ========== PATHS ==========
project_folder = Path ('/scratch/mfafouti/Mommybrain/Slide_tags/Filtering')
in_dir = project_folder / 'NEW_list_merged_filtered'

ad_path = in_dir / "Merged_n=207903_all_metadata_slide_tags.h5ad"
adata = sc.read_h5ad(ad_path)

# Filtering
adata = adata[adata.obs["scDblFinder.class"] == "singlet"].copy()

# Count cells per sample and subclass
counts = adata.obs.groupby(['sample', 'subclass_name']).size().reset_index(name='count')
# Pivot table to make types as rows and samples as columns
pivot = counts.pivot(index='subclass_name', columns='sample', values='count').fillna(0)
types_to_keep = pivot[pivot.columns][(pivot[pivot.columns] >= min_cells).sum(axis=1) >= samples].index.tolist()

print(F"Types to keep: {types_to_keep}")
print("Number of types to keep: ", len(types_to_keep))

# Subset AnnData
adata_filtered = adata[adata.obs['subclass_name'].isin(types_to_keep)].copy()
print(f"Original cells: {adata.n_obs}, Filtered cells: {adata_filtered.n_obs}")

# QC FILTERING
adata_QC_filtered = adata_filtered[adata_filtered.obs['class_bootstrapping_probability']>cutoff_bootstrap].copy()
print(f"Original cells: {adata_filtered.n_obs}, Filtered cells: {adata_QC_filtered.n_obs}")

# DE FILTERING
counts_DE = adata_QC_filtered.obs.groupby(['sample', 'subclass_name']).size().reset_index(name='count')
pivot_DE = counts_DE.pivot(index='subclass_name', columns='sample', values='count').fillna(0)
types_to_keep_DE = pivot_DE[pivot_DE.columns][(pivot_DE[pivot_DE.columns] >= min_cells_DE).all(axis=1)].index.tolist()

print(F"Types to keep: {types_to_keep_DE}")
print("Number of types to keep: ", len(types_to_keep_DE))

adata_filtered_DE = adata_QC_filtered[adata_QC_filtered.obs['subclass_name'].isin(types_to_keep_DE)].copy()
print(f"Original cells: {adata_QC_filtered.n_obs}, Filtered cells: {adata_filtered_DE.n_obs}")

ad_filtered_path = in_dir / f"DE_merged_filtered_{adata_filtered_DE.n_obs}_mincells_{min_cells}_in_{samples}_samples_slide_tags.h5ad"
adata_filtered_DE.write(ad_filtered_path)
print (f'saved object in {ad_filtered_path}')

# # run dim reduction
# print('starting computations')
# sc.pp.highly_variable_genes(adata_QC_filtered, flavor='seurat_v3', n_top_genes=2000)
# sc.pp.normalize_total(adata_QC_filtered, target_sum=1e4)
# sc.pp.log1p(adata_QC_filtered)

# print('done with log normalization')
# sc.pp.scale(adata_QC_filtered, max_value=10)
# sc.tl.pca(adata_QC_filtered, svd_solver='arpack')
# sc.pp.neighbors(adata_QC_filtered, n_neighbors=15, n_pcs=30)
# sc.tl.umap(adata_QC_filtered)
# print('done with umap!')
# adata_QC_filtered.write(os.path.join(in_dir, f"bootsprap_{cutoff_bootstrap}_umap_{adata_QC_filtered.n_obs}_mincells_{min_cells}_in_{samples}_samples_slide_tags.h5ad"))