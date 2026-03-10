import pandas as pd
import numpy as np
from scipy.sparse import issparse
from anndata import AnnData
import matplotlib.pyplot as plt

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

# def collapse_by_gene_symbol(adata, gene_symbol_col="gene_symbols", aggfunc="sum", use_raw=True):
#     """
#     Collapse an AnnData object by gene symbols using raw counts if available.

#     Parameters:
#     -----------
#     adata : AnnData
#         The AnnData object containing gene expression data and gene annotations.
#     gene_symbol_col : str
#         The column in `adata.var` or `adata.raw.var` containing the gene symbols.
#     aggfunc : str or callable
#         Aggregation function to apply for duplicated gene symbols (e.g., 'sum', 'mean').
#     use_raw : bool
#         Whether to use raw counts (adata.raw) if available. Defaults to True.

#     Returns:
#     --------
#     collapsed_adata : AnnData
#         The collapsed AnnData object with unique gene symbols as `var_names` and aggregated expression data.
#     """

#     if use_raw and adata.raw is not None:
#         var_df = adata.raw.var
#         X_matrix = adata.raw.X
#     else:
#         var_df = adata.var
#         X_matrix = adata.X

#     if gene_symbol_col not in var_df.columns:
#         raise ValueError(f"{gene_symbol_col} not found in var annotations")

#     gene_symbols = var_df[gene_symbol_col]

#     # Handle categorical and missing gene symbols in var_df
#     if pd.api.types.is_categorical_dtype(gene_symbols):
#         if "Unknown" not in gene_symbols.cat.categories:
#             gene_symbols = gene_symbols.cat.add_categories(["Unknown"])
#     gene_symbols = gene_symbols.fillna("Unknown")

#     # Filter out unknown gene symbols
#     valid_idx = gene_symbols != "Unknown"
#     X_matrix = X_matrix[:, valid_idx]
#     gene_symbols = gene_symbols[valid_idx]

#     # Convert sparse matrix to dense if needed
#     if issparse(X_matrix):
#         X_matrix = X_matrix.toarray()

#     # Build DataFrame with cells as rows, genes as columns
#     df = pd.DataFrame(X_matrix, index=adata.obs_names, columns=gene_symbols.values)

#     # Collapse genes by gene symbol
#     collapsed_df = df.groupby(axis=1, level=0).aggregate(aggfunc)

#     # Create new AnnData object
#     collapsed_adata = AnnData(
#         X=collapsed_df.values,
#         obs=adata.obs.copy(),
#         var=pd.DataFrame(index=collapsed_df.columns)
#     )
#     collapsed_adata.var_names = collapsed_df.columns

#     return collapsed_adata

# ========== FUNCTION TO PLOT VOLCANO ==========
def plot_volcano(result_dict, title, save=False):
    # Extract DE results
    genes = result_dict['names'][0]
    pvals_adj = result_dict['pvals_adj'][0]
    logfc = result_dict['logfoldchanges'][0]

    # Create DataFrame
    df = pd.DataFrame({
        'gene': genes,
        'logFC': logfc,
        'p_adj': pvals_adj
    })
    df['-log10(p_adj)'] = -np.log10(df['p_adj'])

    # Volcano plot
    plt.figure(figsize=(8, 6))
    plt.scatter(df['logFC'], df['-log10(p_adj)'], s=10, alpha=0.7)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-value')
    plt.title(f'Volcano Plot - {title}')

    # Highlight significant genes
    sig = (df['p_adj'] < 0.05) & (np.abs(df['logFC']) > 1)
    plt.scatter(df.loc[sig, 'logFC'], df.loc[sig, '-log10(p_adj)'],
                color='red', s=10, label='Significant')

    plt.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(-1, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(1, color='gray', linestyle='--', linewidth=0.5)

    plt.legend()
    plt.tight_layout()

    if save:
        plt.savefig(f"volcano_{title.replace(' ', '_')}.png", dpi=300)
