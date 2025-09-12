"""
Title: Extract edgeR input from combined anndata object
Description: Pseudobulk for RCTD "singlets" data by summing raw counts for all samples in each cell type from a single AnnData file
Author:   Maria Eleni Fafouti 
Date: 17-07-2025
"""
import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from collections import defaultdict
import re

# === Input & output paths ===
current_path = '/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR'
input_file = os.path.join(current_path,  "cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad")
output_dir = os.path.join(current_path,"out","pseudobulk_outputs")
os.makedirs(output_dir, exist_ok=True)

# # === Functions ===
# def clean_symbol(name):
#     parts = name.split("-")
#     if parts[0] == ".":   # case where it's just a dot
#         return parts[-1]  # keep after the dash
#     else:
#         return parts[0]
    
# === Load the combined AnnData object ===
adata = sc.read_h5ad(input_file)

# # =========== FIXING THE VAR SECTION ===========
# # Remove the extra gene_id column (keep only the index)
# adata.var = adata.var.drop(columns=["gene_id"])

# # cleaning the gene_symbol column 
# adata.var["gene_symbol"] = adata.var["name"].apply(clean_symbol)

# adata.write('cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad')

# === Keeing only shared celltypes ===
celltype_col = "RCTD_first_type_mouse"  # change if needed

# Get cell types per sample
cts_per_sample = adata.obs.groupby("sample")[celltype_col].unique()

# Find intersection across all samples
common_cts = set.intersection(*map(set, cts_per_sample))
print("Cell types in all samples:", common_cts)

# Subset adata to keep common ones
adata = adata[adata.obs[celltype_col].isin(common_cts)].copy()

adata.var["gene_symbol"] = adata.var["gene_symbol"].astype(str)
adata.var_names = adata.var["gene_symbol"]
adata.var_names_make_unique()

print(adata)

# === Check for required columns ===
if 'RCTD_first_type_mouse' not in adata.obs:
    raise ValueError("'RCTD_first_type_mouse' column not found in adata.obs")
if 'sample' not in adata.obs:
    raise ValueError("'sample' column not found in adata.obs")

# === Collect all genes ===
all_genes = adata.var_names.tolist()

# === Dictionary to hold counts for each cell type ===
pseudobulk_by_celltype = defaultdict(list)  # key = cell type, value = list of (sample_id, summed counts)

# === For each cell type, group by sample and sum counts ===
for celltype in adata.obs['RCTD_first_type_mouse'].unique():
    adata_sub = adata[adata.obs['RCTD_first_type_mouse'] == celltype]
    print(f'Processing celltype: {celltype}...')
    for sample_id in adata_sub.obs['sample'].unique():
        adata_sample = adata_sub[adata_sub.obs['sample'] == sample_id]
        X = adata_sample.X
        summed = X.sum(axis=0)
        if sp.issparse(summed):
            summed = np.asarray(summed).flatten()
        else:
            summed = np.array(summed).flatten()
        pseudobulk_by_celltype[celltype].append((sample_id, summed))
        print(f'  Added sample: {sample_id}')

print('Finished with all cell types!')

# === Write one file per cell type ===
for celltype, sample_sums in pseudobulk_by_celltype.items():
    sample_ids, count_arrays = zip(*sample_sums)
    pseudobulk_matrix = np.vstack(count_arrays).T  # genes x samples
    pseudobulk_df = pd.DataFrame(pseudobulk_matrix, index=all_genes, columns=sample_ids)
    pseudobulk_df.to_csv(f"{output_dir}/{celltype}_counts.tsv", sep="\t")