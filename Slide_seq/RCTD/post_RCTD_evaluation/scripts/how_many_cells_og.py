"""
Title: Plot number of cells per metadata meta
Description:  Parsimonious script that requires 2 arguments: 1) the adata file 2) the metadata meta
Author:   Maria Eleni Fafouti 
Date: 28-08-2025
"""
import argparse
import matplotlib.pyplot as plt
import scanpy as sc
import os

def main():
    # ---- Parse arguments ----
    parser = argparse.ArgumentParser(description="Plot cell counts per two metadata categories from an AnnData object")
    parser.add_argument("--file", "-f", required=True, help="Path to the .h5ad file")
    parser.add_argument("--meta1", "-m1", required=True, help="First metadata column in .obs")
    parser.add_argument("--meta2", "-m2", required=False, help="Second metadata column in .obs")
    parser.add_argument("--out", "-o", required=True, help="Output directory")
    args = parser.parse_args()

    # ---- Load data ----
    adata = sc.read_h5ad(args.file)
    adata.obs_names_make_unique()

    # --- Count cells per combination of meta1 & meta2 ---
    counts = (
        adata.obs
        .groupby([args.meta1])
        .size()
        .reset_index(name="count")
    )

    # # ---- Plot ----
    # plt.figure(figsize=(10, 6))
    # pivot = counts.pivot(index=args.meta1, columns=args.meta2, values="count").fillna(0)
    # pivot.plot(kind="bar", stacked=True, figsize=(10, 6))
    # plt.xlabel(args.meta1)
    # plt.ylabel("Number of cells")
    # plt.title(f"Cell counts by {args.meta1} and {args.meta2}")
    # plt.legend(title=args.meta2, bbox_to_anchor=(1.05, 1), loc="upper left")

    # # ---- Save plot ----
    # os.makedirs(args.out, exist_ok=True)
    # outfile = os.path.join(args.out, f"{args.meta1}_{args.meta2}_counts.png")
    # plt.tight_layout()
    # plt.savefig(outfile, dpi=300)
    # print(f"Plot saved to {outfile}")

    # ---- Save counts to CSV ----
    csv_file = os.path.join(args.out, f"{args.meta1}_counts.csv")
    counts.to_csv(csv_file, index=False)
    print(f"Counts saved to {csv_file}")

if __name__ == "__main__":
    main()

    # # --- Find celltypes that have >=50 cells in ALL samples ---
    # valid_cts = (
    #     counts[counts["count"] >= 10]        # keep only counts >= 50
    #     .groupby(celltype_col)[sample_col]   # check per celltype
    #     .nunique()
    # )
    # # keep only celltypes present in ALL samples with >=50 cells
    # n_samples = adata.obs[sample_col].nunique()
    # valid_cts = valid_cts[valid_cts == n_samples].index

    # print("Cell types with ≥10 cells per sample:", list(valid_cts))

    # # --- Subset the AnnData object ---
    # adata = adata[adata.obs[celltype_col].isin(valid_cts)].copy()