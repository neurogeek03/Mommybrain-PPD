import pandas as pd
import numpy as np
from scipy.sparse import issparse
import os
import anndata as ad
from anndata import AnnData
import scanpy as sc

def filter_1to1_mouse_orthologs(adata, mouse_id_col="mouse_gene_stable_ID"):
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
    n_unknown = (adata.var[mouse_id_col] == "Unknown").sum()
    print(f"❓ {n_unknown} out of {n_total} genes have mouse_gene_id = 'Unknown' ({n_unknown / n_total:.1%})")

    var = adata.var.copy()
    var = var[var[mouse_id_col] != "Unknown"]

    gene_counts = var[mouse_id_col].value_counts()
    n_total_mouse_ids = gene_counts.shape[0]
    n_collapsed = (gene_counts > 1).sum()
    print(f"🔎 Mouse gene IDs total: {n_total_mouse_ids}")
    print(f"🔁 Mouse gene IDs with >1 mapped rat gene: {n_collapsed}")

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