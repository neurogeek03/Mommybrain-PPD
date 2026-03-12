"""
Title: DEA Summary Barplot from LIANA dea_result.csv
Description: Horizontal barplot of up/downregulated gene counts per cell type,
             adapted from 04_plot_edgeR.py to work with LIANA DEA output format.
Author: Maria Eleni Fafouti

Usage
-----
  python plot_dea_barplot.py --input /path/to/dea_result.csv --output barplot.png
  python plot_dea_barplot.py --input /path/to/dea_result.csv --output barplot.png \\
      --group-a PD_CORT --group-b PD_OIL \\
      --logfc-thresh 1.0 --fdr-thresh 0.05

All arguments
-------------
  --input       Path to dea_result.csv (required)
  --output      Output PNG path (required)
  --group-a     Baseline group name, for legend label (default: GroupA)
  --group-b     Comparison group name, for legend label (default: GroupB)
  --logfc-thresh  log2FC threshold (default: 1.0)
  --fdr-thresh    Adjusted p-value threshold (default: 0.05)
  --min-genes     Min DE genes to include a cell type (default: 1)
  --sort-by       Sort by: total | higher | lower (default: total)
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_dea_barplot(input_path, output_path, group_a="GroupA", group_b="GroupB",
                    logfc_thresh=1.0, fdr_thresh=0.05, min_genes=1,
                    sort_by="total", figsize=(12, 9)):

    df = pd.read_csv(input_path, index_col=0)

    # Rename columns to a standard internal format
    df = df.rename(columns={
        "log2FoldChange": "logFC",
        "padj":           "FDR",
    })

    summary = []
    for celltype, group_df in df.groupby("subclass_name"):
        logfc = group_df["logFC"].astype(float).values
        fdr   = group_df["FDR"].fillna(1.0).astype(float).values

        higher = int(np.sum((logfc >  logfc_thresh) & (fdr < fdr_thresh)))
        lower  = int(np.sum((logfc < -logfc_thresh) & (fdr < fdr_thresh)))

        if (higher + lower) >= min_genes:
            summary.append({"subclass": celltype, "higher": higher, "lower": lower})

    if not summary:
        print("No cell types met the threshold criteria.")
        return

    df_summary = pd.DataFrame(summary)
    df_summary["total"] = df_summary["higher"] + df_summary["lower"]
    df_summary = df_summary.sort_values(sort_by)

    label_higher = f"Higher in {group_b}"
    label_lower  = f"Lower in {group_b}"

    fig, ax = plt.subplots(figsize=figsize)

    bar1 = ax.barh(df_summary["subclass"], df_summary["lower"],
                   label=label_lower, color="blue")
    bar2 = ax.barh(df_summary["subclass"], df_summary["higher"],
                   left=df_summary["lower"],
                   label=label_higher, color="red")

    ax.set_xlabel("Number of significant genes")
    ax.set_title(
        f"DE genes per cell type (|log2FC| > {logfc_thresh}, FDR < {fdr_thresh})\n"
        f"{group_b} vs {group_a} (baseline)",
        fontsize=12
    )
    ax.legend(handles=[bar1, bar2])
    plt.tight_layout()

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Horizontal DEA summary barplot from LIANA dea_result.csv."
    )
    parser.add_argument("--input",  "-i", required=True, help="Path to dea_result.csv")
    parser.add_argument("--output", "-o", required=True, help="Output PNG path")
    parser.add_argument("--group-a", default="GroupA", help="Baseline group name (default: GroupA)")
    parser.add_argument("--group-b", default="GroupB", help="Comparison group name (default: GroupB)")
    parser.add_argument("--logfc-thresh", type=float, default=0.1,  help="log2FC threshold (default: 1.0)")
    parser.add_argument("--fdr-thresh",   type=float, default=0.1, help="FDR threshold (default: 0.05)")
    parser.add_argument("--min-genes",    type=int,   default=1,    help="Min DE genes to include a cell type (default: 1)")
    parser.add_argument("--sort-by", choices=["total", "higher", "lower"],
                        default="total", help="Sort bars by this column (default: total)")
    return parser.parse_args()


def main():
    args = parse_args()
    plot_dea_barplot(
        input_path=args.input,
        output_path=args.output,
        group_a=args.group_a,
        group_b=args.group_b,
        logfc_thresh=args.logfc_thresh,
        fdr_thresh=args.fdr_thresh,
        min_genes=args.min_genes,
        sort_by=args.sort_by,
    )


if __name__ == "__main__":
    main()
