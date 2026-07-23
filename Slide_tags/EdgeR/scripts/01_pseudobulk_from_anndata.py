"""
Title: Extract edgeR input from anndata object
Description:  Pseudobulk snRNA-seq data by summing raw counts for all samples in each cell type.
              All genes are passed through; gene-level filtering is delegated to EdgeR (filterByExpr).
Author:   Maria Eleni Fafouti
Date: 09-03-2026
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import argparse
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from my_functions import collapse_by_gene_symbol
from pathlib import Path

# ========== ARGS ==========
parser = argparse.ArgumentParser(description="Pseudobulk scRNA-seq data from AnnData.")
parser.add_argument("--h5ad", type=Path, required=True,
                    help="Path to the input .h5ad file.")
parser.add_argument("--output-dir", type=Path, required=True,
                    help="Directory where per-cell-type *_counts.tsv files will be saved.")
parser.add_argument("--celltype-col", type=str, default="subclass_name",
                    help="obs column to use as cell type label (default: subclass_name).")
parser.add_argument("--gene-symbol-col", type=str, default="gene_symbols",
                    help="var column containing gene symbols (default: gene_symbols).")
parser.add_argument("--score-col", type=str, default=None,
                    help="obs column containing RCTD singlet score (optional).")
parser.add_argument("--spot-class-col", type=str, default=None,
                    help="obs column containing RCTD spot class (optional).")
parser.add_argument("--neuron-score-thresh", type=float, default=None,
                    help="Min singlet score for neurons (Glut/Gaba/IMN cell types). Requires --score-col.")
parser.add_argument("--non-neuron-score-thresh", type=float, default=None,
                    help="Min singlet score for non-neurons (NN cell types). Requires --score-col.")
parser.add_argument("--singlet-only-non-neurons", action="store_true",
                    help="Keep only spot_class==singlet for non-neurons. Requires --spot-class-col.")
args = parser.parse_args()

ad_file = args.h5ad
out_dir = args.output_dir
celltype_col = args.celltype_col
gene_symbol_col = args.gene_symbol_col
score_col = args.score_col
spot_class_col = args.spot_class_col
neuron_score_thresh = args.neuron_score_thresh
non_neuron_score_thresh = args.non_neuron_score_thresh
singlet_only_non_neurons = args.singlet_only_non_neurons

out_dir.mkdir(exist_ok=True, parents=True)

# ========== READING DATA ==========
adata = sc.read_h5ad(ad_file)
print(adata)

adata.raw = adata  # stores current X as raw
adata.X = adata.raw.X.copy()

# # ========== Removing more cells with mt counts ==========
# adata = adata[adata.obs['pct_counts_mt'] < pct_mt_counts].copy()
# print(f'adata object is now subset to {adata.n_obs}')

# ========== RCTD QUALITY FILTERS (optional) ==========
# Non-neurons: cell types ending in _NN; everything else is neuron
if neuron_score_thresh is not None or non_neuron_score_thresh is not None or singlet_only_non_neurons:
    ct = adata.obs[celltype_col].astype(str)
    is_non_neuron = ct.str.endswith('_NN')
    is_neuron = ~is_non_neuron

    keep = pd.Series(True, index=adata.obs_names)

    if neuron_score_thresh is not None:
        if score_col is None:
            raise ValueError("--score-col is required when --neuron-score-thresh is set")
        keep &= ~is_neuron | (adata.obs[score_col] >= neuron_score_thresh)
        print(f"Neuron score filter (>= {neuron_score_thresh}): {keep[is_neuron].sum()} / {is_neuron.sum()} neurons kept")

    if non_neuron_score_thresh is not None:
        if score_col is None:
            raise ValueError("--score-col is required when --non-neuron-score-thresh is set")
        keep &= ~is_non_neuron | (adata.obs[score_col] >= non_neuron_score_thresh)
        print(f"Non-neuron score filter (>= {non_neuron_score_thresh}): {keep[is_non_neuron].sum()} / {is_non_neuron.sum()} non-neurons kept")

    if singlet_only_non_neurons:
        if spot_class_col is None:
            raise ValueError("--spot-class-col is required when --singlet-only-non-neurons is set")
        keep &= ~is_non_neuron | (adata.obs[spot_class_col] == "singlet")
        print(f"Non-neuron singlet filter: {keep[is_non_neuron].sum()} / {is_non_neuron.sum()} non-neurons kept")

    adata = adata[keep].copy()
    print(f"After RCTD filters: {adata.n_obs} cells remaining")

# ========== COLLAPSE ADATA OBJ BY GENE SYMBOL ==========
collapsed_adata = collapse_by_gene_symbol(adata, gene_symbol_col=gene_symbol_col, aggfunc="sum")
print(collapsed_adata)

adata = collapsed_adata

# ========== LOOPING OVER ALL CELL TYPES ==========
for celltype in adata.obs[celltype_col].unique():
    print(f"Processing cell type: {celltype}")
    
    # Subset AnnData to just this cell type
    adata_sub = adata[adata.obs[celltype_col] == celltype]

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