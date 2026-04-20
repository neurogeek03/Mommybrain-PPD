"""
Data loading and extraction utilities for spatial gene expression viewer.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse


def load_adata(path):
    """Load AnnData object from h5ad file."""
    return sc.read_h5ad(path)


def normalize(adata):
    """Return a normalized copy of adata (normalize_total + log1p)."""
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def subset_samples(adata, samples, sample_col="sample"):
    """Return a copy of adata filtered to the specified samples."""
    mask = adata.obs[sample_col].isin(samples)
    return adata[mask].copy()


def validate_genes(adata, genes):
    """
    Check that each gene in `genes` exists in adata.var["gene_symbol"].
    Raises ValueError listing all missing genes.
    """
    available = set(adata.var["gene_symbol"].values)
    missing = [g for g in genes if g not in available]
    if missing:
        raise ValueError(
            f"The following genes were not found in adata.var['gene_symbol']: {missing}"
        )


def get_expression_per_sample(adata, genes, samples, sample_col="sample"):
    """
    Returns dict {gene: {sample: [float list]}}.
    Assumes adata has already been normalized.
    Gene lookup is by gene_symbol column in adata.var.
    """
    # Build a mapping from gene_symbol -> column index
    gene_symbol_arr = adata.var["gene_symbol"].values
    symbol_to_idx = {sym: i for i, sym in enumerate(gene_symbol_arr)}

    result = {}
    for gene in genes:
        col_idx = symbol_to_idx[gene]
        result[gene] = {}
        for sample in samples:
            mask = (adata.obs[sample_col] == sample).values
            # Extract the column from X for the masked rows
            X_sub = adata.X[mask, :]
            if scipy.sparse.issparse(X_sub):
                col_vals = np.asarray(X_sub[:, col_idx].todense()).flatten()
            else:
                col_vals = np.asarray(X_sub[:, col_idx]).flatten()
            result[gene][sample] = col_vals.tolist()

    return result


def get_spatial_per_sample(adata, samples, sample_col="sample", spatial_key="X_spatial"):
    """
    Returns {sample: {"x": [list], "y": [list]}}.
    Spatial coordinates are taken from adata.obsm[spatial_key].
    """
    coords = adata.obsm[spatial_key]
    result = {}
    for sample in samples:
        mask = adata.obs[sample_col] == sample
        mask_arr = np.asarray(mask)
        sub_coords = coords[mask_arr, :]
        result[sample] = {
            "x": sub_coords[:, 0].tolist(),
            "y": sub_coords[:, 1].tolist(),
        }
    return result


def get_celltype_per_sample(adata, samples, sample_col="sample", celltype_col="RCTD_first_type_rat"):
    """
    Returns {sample: [str list]} of cell type labels per sample.
    """
    result = {}
    for sample in samples:
        mask = adata.obs[sample_col] == sample
        result[sample] = adata.obs.loc[mask, celltype_col].astype(str).tolist()
    return result


def get_celltypes(adata, celltype_col="RCTD_first_type_rat"):
    """Return sorted unique list of cell types."""
    return sorted(adata.obs[celltype_col].astype(str).unique().tolist())


def load_celltype_colors(csv_path, term_set="CCN20230722_SUBC", fallback="#808080"):
    """
    Load cell type → hex color mapping from the cluster annotation CSV.

    Filters to the specified term_set level (default: subclass).
    Cleans names the same way as plot_umap.py: replaces spaces, slashes, hyphens with '_'.

    Returns dict {cell_type_label: hex_color}.
    Unmatched cell types will receive `fallback` color at lookup time.
    """
    df = pd.read_csv(csv_path, usecols=["name", "cluster_annotation_term_set_label", "color_hex_triplet"])
    df = df[df["cluster_annotation_term_set_label"] == term_set].copy()
    df["name"] = df["name"].str.replace(r"[ /-]", "_", regex=True)
    return dict(zip(df["name"], df["color_hex_triplet"]))


def load_gene_list(csv_path, adata):
    """
    Load gene symbols from a CSV with a 'gene' column.
    Filters out Ensembl IDs and validates all remaining symbols against adata.
    Returns a sorted list of valid gene symbols.
    """
    df = pd.read_csv(csv_path, usecols=["gene"])
    syms = df["gene"].dropna().astype(str)
    proper = syms[~syms.str.match(r"^ENS[A-Z]+\d+$")].unique().tolist()
    available = set(adata.var["gene_symbol"].values)
    missing = [g for g in proper if g not in available]
    if missing:
        print(f"  Warning: {len(missing)} genes not found in adata, skipping: {missing}")
    return sorted(g for g in proper if g in available)


def get_all_gene_symbols(adata):
    """
    Return sorted list of proper gene symbols, excluding entries that are
    Ensembl IDs (i.e. no gene name was assigned in the reference).
    """
    syms = adata.var["gene_symbol"].dropna().astype(str)
    proper = syms[~syms.str.match(r"^ENS[A-Z]+\d+$")]
    n_proper = proper.nunique()
    n_total = len(adata.var)
    print(f"  {n_proper} genes with proper symbols out of {n_total} total")
    return sorted(proper.unique().tolist())