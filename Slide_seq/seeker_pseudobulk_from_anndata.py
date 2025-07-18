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

# === Input & output paths ===
input_file = "all_slide_seq_singlets_15.h5ad"
output_dir = "pseudobulk_outputs"
os.makedirs(output_dir, exist_ok=True)

# === Load the combined AnnData object ===
adata = sc.read_h5ad(input_file)

# === Check for required columns ===
if 'RCTD_second_type' not in adata.obs:
    raise ValueError("'RCTD_second_type' column not found in adata.obs")
if 'sample' not in adata.obs:
    raise ValueError("'sample' column not found in adata.obs")

# === Collect all genes ===
all_genes = adata.var_names.tolist()

# === Dictionary to hold counts for each cell type ===
pseudobulk_by_celltype = defaultdict(list)  # key = cell type, value = list of (sample_id, summed counts)

# === For each cell type, group by sample and sum counts ===
for celltype in adata.obs['RCTD_second_type'].unique():
    adata_sub = adata[adata.obs['RCTD_second_type'] == celltype]
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