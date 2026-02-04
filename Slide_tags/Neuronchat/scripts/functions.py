import pandas as pd
import numpy as np
from scipy.sparse import issparse
import os
import anndata as ad
from anndata import AnnData
import scanpy as sc
from pathlib import Path

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

def find_file_by_parts(directory, prefix, suffix):
  """
  Finds files in a directory that start with a prefix and end with a suffix.

  Args:
    directory (str): The directory to search in.
    prefix (str): The beginning of the filename.
    suffix (str): The end of the filename.

  Returns:
    A list of Path objects for the files that match.
  """
  search_path = Path(directory)
  # The '*' is a wildcard that matches any number of characters.
  glob_pattern = f"{prefix}*{suffix}"
  
  # rglob searches recursively. Use glob() for a non-recursive search.
  found_files = list(search_path.rglob(glob_pattern))
  
  return found_files

def collapse_by_gene_symbol(adata, gene_symbol_col="gene_symbols", aggfunc="sum"):
    """
    Collapse an AnnData object by gene symbols.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object containing gene expression data and gene annotations.
    gene_symbol_col : str
        The column in `adata.var` containing the gene symbols.
    aggfunc : str or callable
        Aggregation function to apply for duplicated gene symbols (e.g., 'sum', 'mean').

    Returns:
    --------
    collapsed_adata : AnnData
        The collapsed AnnData object with unique gene symbols as `var_names` and aggregated expression data.
    """
    if gene_symbol_col not in adata.var.columns:
        raise ValueError(f"{gene_symbol_col} not found in adata.var")

    # Define gene symbols
    gene_symbols = adata.var[gene_symbol_col]

    # Extract and clean gene symbols (fill NaN with "Unknown" if necessary)
    if pd.api.types.is_categorical_dtype(adata.var[gene_symbol_col]):  # Checks if it's a Categorical
        if "Unknown" not in adata.var[gene_symbol_col].cat.categories:  # Adds "Unknown" as a valid category if it's not already there
            adata.var[gene_symbol_col] = adata.var[gene_symbol_col].cat.add_categories(["Unknown"])

    adata.var[gene_symbol_col] = adata.var[gene_symbol_col].fillna("Unknown")

    # Filter out "Unknown" genes
    valid_idx = adata.var[gene_symbol_col] != "Unknown"
    adata = adata[:, valid_idx].copy()
    gene_symbols = adata.var[gene_symbol_col]

    # Matrix extraction
    X = adata.X.toarray() if issparse(adata.X) else adata.X
    df = pd.DataFrame(X, index=adata.obs_names, columns=gene_symbols.values)

    # Collapse genes by gene symbol (i.e., group columns)
    collapsed_df = df.groupby(axis=1, level=0).aggregate(aggfunc)  # collapse columns

    print(f"Collapsed df shape: {collapsed_df.shape} (should be cells x new_genes)")

    # Now:
    # - collapsed_df rows = cells = obs
    # - collapsed_df columns = new genes = var

    # Create a new AnnData object
    collapsed_adata = AnnData(
        X=collapsed_df.values,
        obs=adata.obs.copy(),  # Cells
        var=pd.DataFrame(index=collapsed_df.columns)  # Genes
    )

    # Update var_names to the collapsed gene symbols
    collapsed_adata.var_names = collapsed_df.columns

    return collapsed_adata