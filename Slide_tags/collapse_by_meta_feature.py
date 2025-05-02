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

# ========== DEFINING ESSENTIAL PATHS ==========
project_path = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags"
adata_dir = os.path.join(project_path, "Post_bender")

def collapse_by_mouse_gene_id(adata, mouse_id_col="mouse_gene_id", aggfunc="sum"):
    """
    Collapse an AnnData object by mouse gene ID.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object containing gene expression data and gene annotations.
    mouse_id_col : str
        The column in `adata.var` containing the mouse gene IDs.
    aggfunc : str or callable
        Aggregation function to apply for duplicated mouse gene IDs (e.g., 'sum', 'mean').

    Returns:
    --------
    collapsed_df : pd.DataFrame
        A DataFrame with mouse_gene_ids as row indices and cells as columns.
    """
    if mouse_id_col not in adata.var.columns:
        raise ValueError(f"{mouse_id_col} not found in adata.var")

    # Define mouse ids 
    mouse_ids = adata.var[mouse_id_col]

    # Extract and clean mouse gene IDs
    # mouse_ids = adata.var[mouse_id_col].fillna("Unknown")
    
    if pd.api.types.is_categorical_dtype(adata.var[mouse_id_col]): # Checks if it's a Categorical
        if "Unknown" not in adata.var[mouse_id_col].cat.categories: # Adds "Unknown" as a valid category if it's not already there
            adata.var[mouse_id_col] = adata.var[mouse_id_col].cat.add_categories(["Unknown"]) # Then fills the NaN values
    
    adata.var[mouse_id_col] = adata.var[mouse_id_col].fillna("Unknown")

    # Filter out "Unknown" genes
    valid_idx = adata.var[mouse_id_col] != "Unknown"
    adata = adata[:, valid_idx].copy()
    mouse_ids = adata.var[mouse_id_col]

    # Matrix extraction
    X = adata.X.toarray() if issparse(adata.X) else adata.X
    df = pd.DataFrame(X, index=adata.obs_names, columns=mouse_ids.values)

    # Collapse genes by mouse_gene_stable_ID (i.e. group columns)
    collapsed_df = df.groupby(axis=1, level=0).aggregate(aggfunc)  # collapse columns

    print(f"Collapsed df shape: {collapsed_df.shape} (should be cells x new_genes)")

    # Now:
    # - collapsed_df rows = cells = obs
    # - collapsed_df columns = new genes = var

    collapsed_adata = AnnData(
        X=collapsed_df.values,
        obs=adata.obs.copy(),  # Cells
        var=pd.DataFrame(index=collapsed_df.columns)  # Genes
    )

    return collapsed_adata


# ========== DEFINING SAMPLE AND TREATMENT LISTS ==========
# sample_list = ["BC13"]
sample_list = ["BC13", "BC14", "BC15", "BC28", "BC3", "BC9"]

for sample in sample_list:
    print(f'Sample {sample} is being processed')
    adata_file = os.path.join(adata_dir, f"{sample}", f"converted_ann_data_{sample}.h5ad")

    # Reading adata file
    adata = sc.read_h5ad(adata_file)
    print(adata.var.head())
    print(adata.var.shape) 
    print(adata.obs.shape) 
    print(adata.obs.head())

    collapsed = collapse_by_mouse_gene_id(adata, mouse_id_col="mouse_gene_stable_ID", aggfunc="sum")
    pd.DataFrame(collapsed.X[:5, :5], index=collapsed.obs_names[:5], columns=collapsed.var_names[:5])


    output_file = os.path.join(adata_dir, f"{sample}", f"collapsed_mouse_id_{sample}.h5ad")
    collapsed.write(output_file)
    print(f"Saved collapsed data to {output_file}")
