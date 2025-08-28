"""
Title: Extract edgeR input from anndata object
Description:  Pseudobulk snRNA-seq data by summing raw counts for all samples in each cell type
Author:   Maria Eleni Fafouti 
Date: 18-06-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import os
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from my_functions import collapse_by_gene_symbol

# ========== DEFINING ESSENTIAL PATHS ==========
project_path = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags'
out_dir = os.path.join(project_path,'EdgeR', 'pseudobulk_output_filtered') 
os.makedirs(out_dir, exist_ok=True)

ad_file = os.path.join(project_path, 'Filtering', 'FILTERED_persample_subclass_doublets_harmony.h5ad')

# ========== READING DATA ==========
adata = sc.read_h5ad(ad_file)
# print(adata)

adata.X = adata.raw.X.copy()

# ========== COLLAPSE ADATA OBJ BY GENE SYMBOL ==========
collapsed_adata = collapse_by_gene_symbol(adata, gene_symbol_col="gene_symbols", aggfunc="mean")
print(collapsed_adata)

adata = collapsed_adata

# ========== LOOPING OVER ALL CELL TYPES ==========
for celltype in adata.obs['MapMyCells_subclass_name'].unique():
    print(f"Processing cell type: {celltype}")
    
    # Subset AnnData to just this cell type
    adata_sub = adata[adata.obs['MapMyCells_subclass_name'] == celltype]

    # === Filter genes: expressed in ≥10% of cells in any condition ===
    conditions = ['CORT', 'OIL']
    keep_genes = set()

    for cond in conditions:
        cond_mask = adata_sub.obs['treatment'] == cond
        X_cond = adata_sub[cond_mask].X

        if sp.issparse(X_cond):
            nonzero_counts = (X_cond > 0).sum(axis=0).A1
        else:
            nonzero_counts = (X_cond > 0).sum(axis=0)

        total_cells = cond_mask.sum()
        frac_expressed = nonzero_counts / total_cells
        gene_mask = frac_expressed >= 0.10
        keep_genes.update(adata_sub.var_names[gene_mask])

    adata_sub = adata_sub[:, list(keep_genes)]
    print(f" {len(keep_genes)} genes retained for cell type '{celltype}' after 10% expression filtering.")

    # Prepare per-sample aggregation
    samples = adata_sub.obs['sample'].unique()
    pseudobulk_data = []
    sample_ids = []

    for sample in samples:
        sample_mask = adata_sub.obs['sample'] == sample
        X = adata_sub[sample_mask].X

        # Sum gene counts across cells for that sample
        summed = X.sum(axis=0)
        if sp.issparse(summed):
            summed = np.asarray(summed).flatten()
        else:
            summed = np.array(summed).flatten()

        pseudobulk_data.append(summed)
        sample_ids.append(sample)
    
    # Create pseudobulk count matrix (genes x samples)
    pseudobulk_matrix = np.vstack(pseudobulk_data).T
    pseudobulk_df = pd.DataFrame(pseudobulk_matrix, index=adata_sub.var_names, columns=sample_ids)

    # === Export both count matrix and sample metadata ===
    safe_celltype = celltype.replace("/", "_").replace(" ", "_")
    pseudobulk_df.to_csv(f"{out_dir}/{safe_celltype}_counts.tsv", sep="\t")

print("✅ All cell types pseudobulked with treatment info.")