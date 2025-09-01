"""
Title: UMI Distribution per Cell
Description: Loads an .h5ad file and plots the distribution of UMIs per cell, using the current matrix, not metadata.
Author: Maria Eleni Fafouti
Date: 01.09.2025
"""

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np

def main():
    # ---- Parse arguments ----
    parser = argparse.ArgumentParser(description="Check UMI distribution per cell from a .h5ad file")
    parser.add_argument("--file", "-f", required=True, help="Path to the .h5ad file")
    parser.add_argument("--outdir", "-o", default="umi_distribution", help="Directory to save plots")
    parser.add_argument("--use_obs", action="store_true", help="Use nCount_RNA from .obs instead of calculating from matrix")
    args = parser.parse_args()

    # ---- Load AnnData ----
    print(f"Loading {args.file} ...")
    adata = sc.read_h5ad(args.file)
    base_name = os.path.splitext(os.path.basename(args.file))[0]

    # ---- Extract UMIs ----
    if args.use_obs:
        if "nCount_RNA" in adata.obs:
            umi_counts = adata.obs["nCount_RNA"].to_numpy()
        elif "total_counts" in adata.obs:
            umi_counts = adata.obs["total_counts"].to_numpy()
        else:
            raise ValueError("No UMI count column found in .obs (expected 'nCount_RNA' or 'total_counts').")
        print("Using UMI counts from .obs")
    else:
        # Calculate directly from the matrix
        umi_counts = adata.X.sum(axis=1)
        # Flatten if sparse or matrix
        if hasattr(umi_counts, "A1"):
            umi_counts = umi_counts.A1
        else:
            umi_counts = np.asarray(umi_counts).ravel()
        print("Using UMI counts calculated from the matrix")

    # ---- Ensure output directory exists ----
    os.makedirs(args.outdir, exist_ok=True)

    # ---- Histogram ----
    plt.figure(figsize=(8,6))
    plt.hist(umi_counts, bins=100, color="skyblue", edgecolor="black")
    plt.axvline(x=100, color='red', linestyle='--', linewidth=2, label='x=100')  # vertical line
    plt.xlabel("UMIs per cell")
    plt.ylabel("Number of cells")
    plt.title(f"Distribution of UMIs per Cell ({base_name})")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f"matrix_{base_name}_umi_distribution_hist.png"), dpi=300)
    plt.close()

    # ---- Boxplot ----
    plt.figure(figsize=(6,4))
    plt.boxplot(umi_counts, vert=False, patch_artist=True)
    plt.xlabel("UMIs per cell")
    plt.title(f"UMI Distribution (Boxplot) - {base_name}")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f"matrix_{base_name}_umi_distribution_box.png"), dpi=300)
    plt.close()

    # ---- Print summary stats ----
    print(f"✅ Plots saved in {args.outdir}")
    print(f"Summary stats:\nMin: {umi_counts.min()}\nMedian: {np.median(umi_counts)}\nMean: {umi_counts.mean()}\nMax: {umi_counts.max()}")

if __name__ == "__main__":
    main()

