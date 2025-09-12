"""
Title:        Custom gene x Cell Matrix from anndata  
Description:  Extracting a gene x cell matrix from an anndata object, using one of the var metadata as the colnames (mouse ID)
Author:       Maria Eleni Fafouti 
Date:         01-05-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import pandas as pd
import numpy as np
from scipy.sparse import issparse
import os
import anndata as ad
from anndata import AnnData
import scanpy as sc

# ========== Paths ==========
project_path = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference_coronal"
adata_dir = os.path.join(project_path, "ABC_regions_mouse_with_rat_genes_h5ads", "full_objects_before_collapse")
out_dir = os.path.join(project_path, "ABC_regions_mouse_with_rat_genes_h5ads")

# ========== Function ==========
def filter_1to1_mouse_orthologs(adata, mouse_id_col="rat_ID"):
    """
    Filter AnnData to keep only genes with 1:1 mouse ortholog mapping.

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object.
    mouse_id_col : str
        Column in adata.var with mouse gene IDs.

    Returns:
    --------
    filtered_adata : AnnData
        New AnnData object filtered to 1:1 mouse ortholog genes.
    """
    if mouse_id_col not in adata.var.columns:
        raise ValueError(f"{mouse_id_col} not found in adata.var")

    n_total = adata.var.shape[0]
    n_unknown = (adata.var[mouse_id_col] == "NaN").sum()
    print(f"❓ {n_unknown} out of {n_total} genes have mouse_gene_id = 'Unknown' ({n_unknown / n_total:.1%})")

    var = adata.var.copy()
    var = var[var[mouse_id_col] != "Unknown"]

    gene_counts = var[mouse_id_col].value_counts()
    n_total_mouse_ids = gene_counts.shape[0]
    n_collapsed = (gene_counts > 1).sum()
    print(f"🔎 {mouse_id_col} gene IDs total: {n_total_mouse_ids}")
    print(f"🔁 {mouse_id_col} gene IDs with >1 mapped mouse gene: {n_collapsed}")

    unique_mouse_genes = gene_counts[gene_counts == 1].index
    is_1to1 = var[mouse_id_col].isin(unique_mouse_genes)

    filtered_genes = var[is_1to1].index
    final_mask = adata.var.index.isin(filtered_genes)

    filtered_adata = adata[:, final_mask].copy()
    filtered_adata.var.set_index(mouse_id_col, inplace=True)

    #Sense checks: 
    print(f"✅ Filtered AnnData shape: {filtered_adata.shape} (cells × genes)")
    print("📋 Head of filtered var:")
    print(filtered_adata.var.head(20))
    
    return filtered_adata


# ========== Running ==========
for fname in os.listdir(adata_dir):
    if not fname.endswith(".h5ad"):
        continue

    fpath = os.path.join(adata_dir, fname)
    print(f"🔄 Processing {fname}...")

    adata = ad.read_h5ad(fpath)

    print(adata.var.head())
    print(adata.var.shape) 
    print(adata.obs.shape) 
    print(adata.obs.head())

    collapsed = filter_1to1_mouse_orthologs(adata, mouse_id_col="rat_ID")
    basename = os.path.basename(fname)
    print(f"The name of this file will start with {basename}")

    output_file = os.path.join(out_dir, f"{basename}_collapsed_rat_genes.h5ad")
    collapsed.write(output_file)
    print(f"Saved collapsed data to {output_file}")
