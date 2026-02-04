from anndata import AnnData
import pandas as pd
import numpy as np
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

def print_filter_debug_info(ctdata: AnnData, cell_group: str, condition_key: str):
    """
    Prints debugging information for the ctdata object before filtering.

    Parameters:
    -----------
    ctdata : anndata.AnnData
        The AnnData object for the current cell group.
    cell_group : str
        The name of the current cell group.
    condition_key : str
        The key for the condition in ctdata.obs (e.g., 'group').
    """
    print(f"\n--- Debugging for cell group: {cell_group} ---")
    print("ctdata object:")
    print(ctdata)
    print("\nctdata.obs:")
    print(ctdata.obs)
    print(f"\nValue counts for '{condition_key}' column:")
    print(ctdata.obs[condition_key].value_counts())
    print("\nData type of the matrix:")
    print(ctdata.X.dtype)

    gene_to_check = '6330411D24Rik'
    if gene_to_check in ctdata.var_names:
        print(f"\nCounts for gene '{gene_to_check}':")
        gene_idx = ctdata.var_names.get_loc(gene_to_check)
        # Handle sparse vs dense for printing
        if hasattr(ctdata.X, 'toarray'):
            print(ctdata.X[:, gene_idx].toarray().flatten())
        else:
            print(ctdata.X[:, gene_idx])
    else:
        print(f"\nGene '{gene_to_check}' not found in this ctdata. Current var_names: {ctdata.var_names[:5].tolist()}...")
    print("--- End of debug block ---")

def preview_X(matrix_like, n_rows=5, n_cols=5, row_labels=None, col_labels=None, title="Matrix Preview"):
    """
    Generates a formatted string preview of a matrix-like object
    (e.g., numpy array, scipy sparse matrix, AnnData.X).

    Parameters:
    -----------
    matrix_like : np.ndarray, sp.spmatrix, or AnnData.X
        The matrix or matrix-like object to preview.
    n_rows : int
        Number of rows to display in the preview.
    n_cols : int
        Number of columns to display in the preview.
    row_labels : list-like, optional
        Labels for the rows (e.g., obs_names).
    col_labels : list-like, optional
        Labels for the columns (e.g., var_names).
    title : str
        Title for the preview.

    Returns:
    --------
    str
        A formatted string representation of the matrix preview.
    """

    if matrix_like is None:
        return f"{title}:\n    (Matrix is None)"

    original_shape = getattr(matrix_like, "shape", (0, 0))
    if original_shape == (0, 0):
        return f"{title}:\n    (Matrix is empty)"

    # Handle sparse matrices
    if issparse(matrix_like):
        matrix_dense = matrix_like.toarray()
    else:
        matrix_dense = matrix_like

    # Get subset for preview
    preview_matrix = matrix_dense[
        :min(n_rows, original_shape[0]),
        :min(n_cols, original_shape[1])
    ]

    # Get labels
    preview_row_labels = (
        row_labels[:min(n_rows, original_shape[0])]
        if row_labels is not None
        else [f"row_{i}" for i in range(preview_matrix.shape[0])]
    )

    preview_col_labels = (
        col_labels[:min(n_cols, original_shape[1])]
        if col_labels is not None
        else [f"col_{i}" for i in range(preview_matrix.shape[1])]
    )

    df_preview = pd.DataFrame(
        preview_matrix,
        index=preview_row_labels,
        columns=preview_col_labels
    )

    header = f"--- {title} (Shape: {original_shape}) ---"
    footer = "----------------------------------"

    return f"{header}\n{df_preview.to_string()}\n{footer}"

