"""
Title: Plotting EdgeR results
Description: Comparing the expression of a pseudobulked cell type across animals
Author: Maria Eleni Fafouti
Date: 23-06-2025

Usage
-----
Two subcommands are available: `barplot` and `volcanos`.

  # Summary barplot of DE gene counts per subclass (default thresholds: |log2FC| > 1, FDR < 0.05)
  python plot_edgeR.py barplot

  # With custom thresholds, horizontal bars, sorted by upregulated count
  python plot_edgeR.py barplot --logfc-thresh 0.5 --fdr-thresh 0.1 --horizontal --sort-by upregulated

  # Individual volcano plots for each subclass (default thresholds: |log2FC| > 0.1, FDR < 0.1)
  python plot_edgeR.py volcanos

  # With custom input/output directories
  python plot_edgeR.py --input-dir /path/to/edger_results volcanos --output-dir /path/to/out

All arguments
-------------
  --input-dir / -i    Directory with *_edgeR_results.tsv files (shared across subcommands)

  barplot:
    --output / -o     Output file path for the barplot PNG
    --logfc-thresh    log2FC threshold (default: 1.0)
    --fdr-thresh      FDR threshold (default: 0.05)
    --min-genes       Min DE genes required to include a subclass (default: 1)
    --sort-by         Sort bars by: total | upregulated | downregulated (default: total)
    --horizontal      Plot horizontal bars instead of vertical

  volcanos:
    --output-dir / -o Directory to save individual volcano PNGs
    --logfc-thresh    log2FC threshold (default: 0.1)
    --fdr-thresh      FDR threshold (default: 0.1)
"""

# ========== IMPORTS ==========
import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ========== PARAMETERS ==========
PROJECT_DIR = Path.cwd().parents[0]
DEFAULT_INPUT_DIR = os.path.join(PROJECT_DIR, "out", "new_march_26", "3pct_mt_pct", "edger_lrt")
DEFAULT_OUTPUT_DIR = os.path.join(PROJECT_DIR, "out", "figures")

# Gene labels for specific subclasses (used by volcanos command)
SUBCLASS_GENES = {
    "327_Oligo_NN": ["Fkbp5", "Sgk1", "Ddit4", "Pdk4", "Btg2"],
    "007_L2_3_IT_CTX_Glut": ["Fkbp5", "Btg2", "Zfp189", "Tgfb2", "Wipf3"],
    "061_STR_D1_Gaba": ["Camk1g", "Actg1", "Kcnh7", "Arhgap26"],
}


# ========== FUNCTIONS ==========
def plot_volcano_edger(res_df, title, ax=None, logfc_thresh=0.1, fdr_thresh=0.1):
    """
    Plots a volcano plot for differential expression analysis results.

    Parameters:
    - res_df: DataFrame with 'logFC', 'PValue', and 'FDR' columns
    - title: str, plot title
    - ax: matplotlib Axes, optional
    - logfc_thresh: float, log2 fold change threshold
    - fdr_thresh: float, FDR threshold
    """
    logfc = res_df["logFC"].astype(float).values
    pval = res_df["PValue"].astype(float).values
    fdr = res_df["FDR"].astype(float).values
    genes = res_df.index.values

    valid = ~np.isnan(logfc) & ~np.isnan(pval) & ~np.isnan(fdr)
    logfc, pval, fdr, genes = logfc[valid], pval[valid], fdr[valid], genes[valid]

    sig = (np.abs(logfc) > logfc_thresh) & (fdr < fdr_thresh)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7))

    ax.scatter(logfc, -np.log10(pval), c="gray", alpha=0.6, label="Not sig.")
    ax.scatter(logfc[sig], -np.log10(pval[sig]), c="#732B8B", alpha=0.8, label=f"FDR < {fdr_thresh}")

    ax.axvline(-logfc_thresh, color="black", linestyle="dotted", lw=1)
    ax.axvline(logfc_thresh, color="black", linestyle="dotted", lw=1)
    ax.axhline(-np.log10(0.05), color="black", linestyle="dashed", lw=1)

    ax.set_title(title, fontsize=14)
    ax.set_xlabel("log2 fold change", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.legend()


def plot_deg_barplot(input_dir, output_path, logfc_thresh=1, fdr_thresh=0.5,
                     min_genes=1, sort_by="total", figsize=(12, 9), horizontal=False):
    """
    Generate a barplot of number of significantly up/downregulated genes per file.

    Parameters:
    - input_dir: str, directory with *_dge_results.tsv files
    - output_path: str, where to save the plot
    - logfc_thresh: float, log2 fold change threshold
    - fdr_thresh: float, FDR cutoff
    - min_genes: int, minimum number of DE genes to include a group
    - sort_by: str, one of "total", "upregulated", "downregulated"
    - figsize: tuple, size of the figure
    - horizontal: bool, if True, plot horizontal bars
    """
    file_paths = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith("_edgeR_results.tsv")
    ]
    summary = []

    for filepath in sorted(file_paths):
        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)
            logfc = df["logFC"].astype(float).values
            fdr = df["FDR"].astype(float).values

            up = np.sum((logfc > logfc_thresh) & (fdr < fdr_thresh))
            down = np.sum((logfc < -logfc_thresh) & (fdr < fdr_thresh))

            if (up + down) >= min_genes:
                title = os.path.basename(filepath).replace("_edgeR_results.tsv", "")
                summary.append({"subclass": title, "upregulated": up, "downregulated": down})

        except Exception as e:
            print(f"Skipping {filepath} due to error: {e}")

    if not summary:
        print("No valid DE files found or no DE genes met threshold.")
        return

    df_summary = pd.DataFrame(summary)
    df_summary["total"] = df_summary["upregulated"] + df_summary["downregulated"]
    df_summary = df_summary.sort_values(sort_by)

    fig, ax = plt.subplots(figsize=figsize)
    if horizontal:
        bar1 = ax.barh(df_summary["subclass"], df_summary["downregulated"],
                       label="Downregulated (OIL > CORT)", color="blue")
        bar2 = ax.barh(df_summary["subclass"], df_summary["upregulated"],
                       left=df_summary["downregulated"],
                       label="Upregulated (CORT > OIL)", color="red")
        ax.set_xlabel("Number of significant genes")
    else:
        x = np.arange(len(df_summary))
        width = 0.35
        bar1 = ax.bar(x - width / 2, df_summary["downregulated"], width,
                      label="Downregulated (OIL > CORT)", color="blue")
        bar2 = ax.bar(x + width / 2, df_summary["upregulated"], width,
                      label="Upregulated (CORT > OIL)", color="red")
        ax.set_xticks(x)
        ax.set_xticklabels(df_summary["subclass"], rotation=90, fontsize=8)
        ax.set_ylabel("Number of significant genes")

    ax.set_title(f"DE genes per subclass (|log2FC| > {logfc_thresh}, FDR < {fdr_thresh})", fontsize=12)
    ax.legend(handles=[bar1, bar2])
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")


def plot_and_save_individual_volcanos(input_dir, output_dir, logfc_thresh=0.1, fdr_thresh=0.1):
    """
    Creates and saves individual volcano plots for each DE result file,
    labeling subclass-specific genes defined in SUBCLASS_GENES.

    Parameters:
    - input_dir: str, directory containing *_dge_results.tsv files
    - output_dir: str, directory where individual plots will be saved
    - logfc_thresh: float, log2 fold change threshold for significance
    - fdr_thresh: float, FDR threshold for significance
    """
    os.makedirs(output_dir, exist_ok=True)

    file_paths = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith("_edgeR_results.tsv")
    ]

    for filepath in sorted(file_paths):
        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)
            df = df.dropna(subset=["logFC", "PValue", "FDR"])

            title = os.path.basename(filepath).replace("_edgeR_results.tsv", "")
            fig, ax = plt.subplots(figsize=(8, 6))

            plot_volcano_edger(df, title=title, ax=ax,
                               logfc_thresh=logfc_thresh, fdr_thresh=fdr_thresh)

            if title in SUBCLASS_GENES:
                for gene in SUBCLASS_GENES[title]:
                    if gene in df.index:
                        row = df.loc[gene]
                        ax.annotate(
                            gene,
                            (row["logFC"], -np.log10(row["PValue"])),
                            fontsize=10,
                            xytext=(4, 4),
                            textcoords="offset points",
                            arrowprops=dict(arrowstyle="-", lw=0.5),
                        )

            plt.tight_layout()
            outpath = os.path.join(output_dir, f"{title}_volcano.png")
            plt.savefig(outpath, dpi=300)
            plt.close()
            print(f"Saved: {outpath}")

        except Exception as e:
            print(f"Skipping {filepath} due to error: {e}")


# ========== CLI ==========
def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot EdgeR differential expression results."
    )
    parser.add_argument(
        "--input-dir", "-i",
        default=DEFAULT_INPUT_DIR,
        help=f"Directory with *_edgeR_results.tsv files (default: {DEFAULT_INPUT_DIR})",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- barplot ---
    bp = subparsers.add_parser("barplot", help="Summary barplot of DE gene counts per subclass.")
    bp.add_argument("--output", "-o", default=os.path.join(DEFAULT_OUTPUT_DIR, "DE_summary_barplot.png"),
                    help="Output file path.")
    bp.add_argument("--logfc-thresh", type=float, default=1, help="log2FC threshold (default: 1.0)")
    bp.add_argument("--fdr-thresh", type=float, default=0.1, help="FDR threshold (default: 0.05)")
    bp.add_argument("--min-genes", type=int, default=1, help="Min DE genes to include a group (default: 1)")
    bp.add_argument("--sort-by", choices=["total", "upregulated", "downregulated"],
                    default="total", help="Sort bars by this column (default: total)")
    bp.add_argument("--horizontal", action="store_true", help="Plot horizontal bars.")

    # --- volcanos ---
    vp = subparsers.add_parser("volcanos", help="Individual volcano plots per subclass.")
    vp.add_argument("--output-dir", "-o", default=DEFAULT_OUTPUT_DIR,
                    help=f"Directory for output plots (default: {DEFAULT_OUTPUT_DIR})")
    vp.add_argument("--logfc-thresh", type=float, default=0.1, help="log2FC threshold (default: 0.1)")
    vp.add_argument("--fdr-thresh", type=float, default=0.1, help="FDR threshold (default: 0.1)")

    return parser.parse_args()


def main():
    args = parse_args()

    if args.command == "barplot":
        plot_deg_barplot(
            input_dir=args.input_dir,
            output_path=args.output,
            logfc_thresh=args.logfc_thresh,
            fdr_thresh=args.fdr_thresh,
            min_genes=args.min_genes,
            sort_by=args.sort_by,
            horizontal=args.horizontal,
        )

    elif args.command == "volcanos":
        plot_and_save_individual_volcanos(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            logfc_thresh=args.logfc_thresh,
            fdr_thresh=args.fdr_thresh,
        )


if __name__ == "__main__":
    main()
