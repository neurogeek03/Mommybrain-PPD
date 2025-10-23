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

pct_mt_counts = 5

# paths
project_folder = Path ('/scratch/mfafouti/Mommybrain/Slide_tags/Filtering')
in_dir = project_folder / 'NEW_list_merged_filtered'
out_dir = project_folder / 'out'
out_dir.mkdir(exist_ok=True)
fig_path = out_dir / 'QC_plot.png'

ad_path = in_dir / "Merged_n=207903_all_metadata_slide_tags.h5ad"
adata = sc.read_h5ad(ad_path)

# ============= DOUBLETS =============
adata = adata[adata.obs["scDblFinder.class"] == "singlet"].copy()

# ============= UMI filtering =============
# sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# ============= MITOCHONDRIA =============
adata.var['mt'] = adata.var['gene_symbol'].str.startswith('Mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p = False)

# sc.pl.violin(
#     adata,
#     ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
#     jitter=0.4,
#     multi_panel=True,
#     save = fig_path
# )

# pd.set_option('display.max_columns', None)
# print(adata.var.head())
# print(adata.obs.head())
# print(adata.var['mt'].unique())
# print(adata.obs['n_genes_by_counts'].unique())
# print(adata.obs['total_counts'].unique())
# print(adata.obs['pct_counts_mt'].unique())
# print(adata.obs['pct_counts_mt'].unique())

sum_mt_cells = (adata.obs['pct_counts_mt'] < pct_mt_counts).sum()
print(f'{sum_mt_cells} cells have more than {pct_mt_counts}% mitochondrial genes')
# subsetting
adata = adata[adata.obs['pct_counts_mt'] < pct_mt_counts].copy()

# ============= PER SUBCLASS FILTERING =============
# Count cells per sample and subclass
counts = adata.obs.groupby(['sample', 'subclass_name']).size().reset_index(name='count')
# Pivot table to make types as rows and samples as columns
pivot = counts.pivot(index='subclass_name', columns='sample', values='count').fillna(0)
types_to_keep = pivot[pivot.columns][(pivot[pivot.columns] >= min_cells).sum(axis=1) >= samples].index.tolist()
print(F"Types to keep: {types_to_keep}")
print("Number of types to keep: ", len(types_to_keep))
# Subset AnnData
adata = adata[adata.obs['subclass_name'].isin(types_to_keep)].copy()
print(f"Original cells: {adata.n_obs}, Filtered cells: {adata.n_obs}")

# ============= QC FILTERING =============
adata = adata[adata.obs['class_bootstrapping_probability']>cutoff_bootstrap].copy()
print(f"Original cells: {adata.n_obs}, Filtered cells: {adata.n_obs}")

# ============= DE FILTERING =============
counts_DE = adata.obs.groupby(['sample', 'subclass_name']).size().reset_index(name='count')
pivot_DE = counts_DE.pivot(index='subclass_name', columns='sample', values='count').fillna(0)
types_to_keep_DE = pivot_DE[pivot_DE.columns][(pivot_DE[pivot_DE.columns] >= min_cells_DE).all(axis=1)].index.tolist()

print(F"Types to keep: {types_to_keep_DE}")
print("Number of types to keep: ", len(types_to_keep_DE))

adata_DE = adata[adata.obs['subclass_name'].isin(types_to_keep_DE)].copy()
print(f"Original cells: {adata.n_obs}, Filtered cells: {adata_DE.n_obs}")

# # ============= DIM REDUCTION & UMAP =============
# print('starting computations')
# sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# print('done with log normalization')
# sc.pp.scale(adata, max_value=10)
# sc.tl.pca(adata, svd_solver='arpack')
# sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
# sc.tl.umap(adata)
# print('done with umap!')

# ============= SAVING =============
adata_to_save = adata_DE.copy()
ad_filtered_path = out_dir / f"DE_after_mt_filter_{adata_to_save.n_obs}_mincells_{min_cells}_in_{samples}_samples_slide_tags.h5ad"
adata_to_save.write(ad_filtered_path)
print (f'saved object in {ad_filtered_path}')


# scratch

# MITOCHONDRIAL & PROTEIN CODING FILTERING
# is_protein_coding = adata.var['mouse_protein_transcript_stable_ID'].notna() & (adata.var['Protein stable ID'] != "")

# # Calculate how many mitochondrial genes are expressed in each cell
# mt_expr = (adata[:, mt_genes].X > 0).sum(axis=1).A1  # .A1 flattens sparse matrices
# num_cells_with_mt = (mt_expr > 0).sum()
# print(f"Number of cells with at least one mitochondrial gene: {num_cells_with_mt}")
# adata.obs['n_mt_genes_expressed'] = mt_expr

# is_not_mitochondrial = ~adata.var['gene_symbol'].str.startswith('Mt-')

# adata_mt_removed = adata[:, is_not_mitochondrial].copy()
# print(f"Original genes: {adata.n_vars}, Filtered genes: {adata_mt_removed.n_vars}")
