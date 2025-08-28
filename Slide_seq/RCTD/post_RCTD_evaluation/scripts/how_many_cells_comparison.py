import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd

def main():
    # ---- Parse arguments ----
    parser = argparse.ArgumentParser(description="Plot cell counts per metadata category from an AnnData object")
    parser.add_argument("--file", "-f", required=True, help="Path to the .h5ad file")
    parser.add_argument("--meta", "-m", required=True, help="Metadata column in .obs to count")
    parser.add_argument("--out", "-o", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ---- Load data ----
    adata = sc.read_h5ad(args.file)
    adata.obs_names_make_unique()

    celltype_col = "RCTD_first_type_mouse"  # adjust if needed
    sample_col = "sample"                   # adjust if needed

    # ---- Count cells for the requested metadata (unfiltered) ----
    counts_all = adata.obs[args.meta].value_counts().reset_index()
    counts_all.columns = [args.meta, "count"]

    # ---- Count cells per sample and celltype for filtering ----
    counts_per_sample = (
        adata.obs
        .groupby([sample_col, celltype_col])
        .size()
        .reset_index(name="count")
    )

    # # --- Find celltypes that have >=10 cells in ALL samples ---
    # valid_cts = (
    #     counts_per_sample[counts_per_sample["count"] >= 10]
    #     .groupby(celltype_col)[sample_col]
    #     .nunique()
    # )
    # n_samples = adata.obs[sample_col].nunique()
    # valid_cts = valid_cts[valid_cts == n_samples].index
    # print("Cell types with ≥10 cells per sample:", list(valid_cts))

    # # --- Subset the AnnData object for filtered celltypes ---
    # adata_filtered = adata[adata.obs[celltype_col].isin(valid_cts)].copy()
    # counts_filtered = adata_filtered.obs[args.meta].value_counts().reset_index()
    # counts_filtered.columns = [args.meta, "count"]

    # ---- Save counts to CSV ----
    csv_file_all = os.path.join(args.out, f"{args.meta}_counts_all.csv")
    csv_file_filtered = os.path.join(args.out, f"{args.meta}_counts_filtered.csv")
    counts_all.to_csv(csv_file_all, index=False)
    counts_filtered.to_csv(csv_file_filtered, index=False)
    print(f"Counts saved to {csv_file_all} and {csv_file_filtered}")

    # ---- Plot ----
    plt.figure(figsize=(6, 8))
    ax = counts_all.set_index(args.meta)["count"].plot(kind="barh", color="lightgray", label="All cells")
    counts_filtered.set_index(args.meta)["count"].plot(kind="barh", color="steelblue", label="≥10 per sample", ax=ax)
    plt.xlabel("Number of cells")
    plt.ylabel(args.meta)
    plt.title(f"Cell counts per {args.meta}")
    plt.legend()

    # ---- Save plot ----
    plot_file = os.path.join(args.out, f"{args.meta}_counts_comparison.png")
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    print(f"Plot saved to {plot_file}")

if __name__ == "__main__":
    main()
