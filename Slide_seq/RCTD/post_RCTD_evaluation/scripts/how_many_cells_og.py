"""
Title: Plot number of cells per metadata meta
Description:  Parsimonious script that requires 2 arguments: 1) the adata file 2) the metadata meta
Author:   Maria Eleni Fafouti 
Date: 23-06-2025
"""
import argparse
import matplotlib.pyplot as plt
import scanpy as sc
import os

def main():
    # ---- Parse arguments ----
    parser = argparse.ArgumentParser(description="Plot cell counts per metadata category from an AnnData object")
    parser.add_argument("--file", "-f", required=True, help="Path to the .h5ad file")
    parser.add_argument("--meta", "-m", required=True, help="Metadata meta in .obs to count")
    parser.add_argument("--out", "-o", required=True, help="Output dir")
    args = parser.parse_args()

    # ---- Load data ----
    adata = sc.read_h5ad(args.file)

    adata.obs_names_make_unique()

    celltype_col = "RCTD_first_type_mouse"  # metadata column
    sample_col = "sample"                   # adjust if needed

    # --- Count cells per sample and celltype ---
    counts = (
        adata.obs
        .groupby([sample_col, celltype_col])
        .size()
        .reset_index(name="count")
    )

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

    # ---- Count cells ----
    counts = adata.obs[args.meta].value_counts()

    # ---- Plot ----
    plt.figure(figsize=(6, 8))
    counts.plot(kind="barh")   # <-- horizontal bar plot instead of bar
    plt.xlabel("Number of cells")
    plt.ylabel(args.meta)    # or your metadata column name
    plt.title(f"Cell counts per {args.meta}")

    # ---- Save ----
    outfile = f"{args.meta}_counts.png"
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"Plot saved to {args.out}")

    # ---- Count cells for the requested metadata ----
    counts_meta = adata.obs[args.meta].value_counts().reset_index()
    counts_meta.columns = [args.meta, "count"]

    # ---- Save counts to CSV ----
    csv_file = os.path.join(args.out, f"{args.meta}_counts.csv")
    counts_meta.to_csv(csv_file, index=False)
    print(f"Counts saved to {csv_file}")

if __name__ == "__main__":
    main()
