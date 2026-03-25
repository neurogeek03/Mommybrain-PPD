"""
Extract cell type counts and sample metadata from h5ad for crumblr.

Outputs:
  <out_dir>/cell_counts.csv    — rows=samples, cols=cell types (subclass name)
  <out_dir>/sample_metadata.csv — rows=samples, col=treatment
"""

import argparse
import os

import pandas as pd
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True, help="Path to h5ad file")
    p.add_argument("--sample_col", default="sample",
                   help="obs column for sample ID")
    p.add_argument("--celltype_col", default="subclass name",
                   help="obs column for cell type")
    p.add_argument("--treatment_col", default="treatment",
                   help="obs column for treatment")
    p.add_argument("--out_dir", required=True, help="Output directory")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading {args.h5ad} ...")
    adata = sc.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

    # --- cell counts: samples x cell types ---
    counts = (
        adata.obs
        .groupby([args.sample_col, args.celltype_col], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    counts.index.name = "sample"
    counts_path = os.path.join(args.out_dir, "cell_counts.csv")
    counts.to_csv(counts_path)
    print(f"Saved cell counts: {counts.shape[0]} samples x {counts.shape[1]} cell types -> {counts_path}")

    # --- sample metadata: one row per sample ---
    metadata = (
        adata.obs[[args.sample_col, args.treatment_col]]
        .drop_duplicates(subset=args.sample_col)
        .set_index(args.sample_col)
        .sort_index()
    )
    meta_path = os.path.join(args.out_dir, "sample_metadata.csv")
    metadata.to_csv(meta_path)
    print(f"Saved sample metadata: {metadata.shape[0]} samples -> {meta_path}")
    print(f"  Treatment counts:\n{metadata[args.treatment_col].value_counts().to_string()}")


if __name__ == "__main__":
    main()