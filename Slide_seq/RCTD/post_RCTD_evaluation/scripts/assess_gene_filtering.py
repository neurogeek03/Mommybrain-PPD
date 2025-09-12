#!/usr/bin/env python3
"""
Assess sparse expression filtering (counts_MIN, MIN_OBS) on a .h5ad file.

Usage:
    python assess_gene_filtering.py --file data.h5ad --counts-min 10 --min-obs 3 --outdir qc_plots
"""

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    # ---- Parse arguments ----
    parser = argparse.ArgumentParser(description="Assess RCTD-style gene filtering on a .h5ad file")
    parser.add_argument("--file", "-f", required=True, help="Path to the .h5ad file")
    parser.add_argument("--counts-min", type=int, default=10, help="Minimum total counts per gene (default: 10)")
    parser.add_argument("--min-obs", type=int, default=3, help="Minimum number of cells expressing a gene (default: 3)")
    parser.add_argument("--outdir", "-o", default="qc_plots", help="Directory to save plots")
    args = parser.parse_args()

    # ---- Load data ----
    print(f"Loading {args.file} ...")
    adata = sc.read_h5ad(args.file)

    counts = adata.X
    if not hasattr(counts, "tocsc"):  # ensure sparse handling
        counts = sc.AnnData(counts).X

    # ---- Gene-level filtering ----
    gene_totals = np.array(counts.sum(axis=0)).flatten()
    gene_obs = np.array((counts > 0).sum(axis=0)).flatten()

    keep_genes = (gene_totals >= args.counts_min) & (gene_obs >= args.min_obs)

    print(f"Genes before filtering: {len(gene_totals)}")
    print(f"Genes after filtering : {np.sum(keep_genes)}")
    print(f"Genes dropped         : {np.sum(~keep_genes)}")

    # ---- Cell-level impact ----
    filtered_counts = counts[:, keep_genes]
    genes_per_cell = np.array((filtered_counts > 0).sum(axis=1)).flatten()

    print(f"Cells: {filtered_counts.shape[0]}")
    print(f"Median genes per cell after filtering: {np.median(genes_per_cell):.1f}")

    # ---- Make output dir ----
    os.makedirs(args.outdir, exist_ok=True)

    # ---- Plots ----
    # Histogram of genes per cell
    plt.figure(figsize=(6,4))
    plt.hist(genes_per_cell, bins=50, color="steelblue", edgecolor="black")
    plt.xlabel("Detected genes per cell (after filtering)")
    plt.ylabel("Number of cells")
    plt.title("Genes per cell distribution")
    plt.axvline(50, color="red", linestyle="--", label="50-gene cutoff")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "genes_per_cell_hist.png"))
    plt.close()

    # Histogram of total counts per gene before vs after
    plt.figure(figsize=(6,4))
    plt.hist(np.log10(gene_totals + 1), bins=50, alpha=0.5, label="Before filtering")
    plt.hist(np.log10(gene_totals[keep_genes] + 1), bins=50, alpha=0.5, label="After filtering")
    plt.xlabel("log10(total counts per gene)")
    plt.ylabel("Number of genes")
    plt.title("Gene filtering effect")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "gene_totals_hist.png"))
    plt.close()

    print(f"Plots saved to {args.outdir}")

if __name__ == "__main__":
    main()
