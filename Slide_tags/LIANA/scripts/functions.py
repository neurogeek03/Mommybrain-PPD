import anndata
import pandas as pd
import adata as ad
import scipy.sparse as sp
from scipy.sparse import issparse

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