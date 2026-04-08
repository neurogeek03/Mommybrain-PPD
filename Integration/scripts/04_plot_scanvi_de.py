"""
04_plot_scanvi_de.py
--------------------
Diverging barplot of scANVI DE results: genes higher in CORT vs OIL per cell type.

FDR is computed per cell type by applying BH correction to (1 - proba_m1) for
upregulated genes and (1 - proba_m2) for downregulated genes, then filtering on
fdr_thresh and |lfc| > lfc_thresh.

Reads per-cell-type CSV files from the DE output directory.
Plotting logic adapted from Slide_tags/EdgeR/scripts/04_plot_edgeR.py.

Usage:
  python 04_plot_scanvi_de.py
  python 04_plot_scanvi_de.py --input-dir out/configB/de/per_cell_type
                              --output out/configB/de/de_barplot.png
                              --fdr-thresh 0.05 --lfc-thresh 0.25 --min-genes 5
"""

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# =============================================================================
# DEFAULT PATHS
# =============================================================================

DEFAULT_INPUT_DIR = "out/configB/de/per_cell_type"
DEFAULT_OUTPUT    = "out/configB/de/de_barplot.png"


# =============================================================================
# PLOT FUNCTION (adapted from plot_deg_barplot in 04_plot_edgeR.py)
# =============================================================================

def plot_scanvi_de_barplot(
    input_dir,
    output_path,
    fdr_thresh=0.05,
    lfc_thresh=0.25,
    min_genes=1,
    sort_by="total",
    figsize=(8, 20),
    group1="CORT",
    group2="OIL",
):
    """
    Diverging barplot of scANVI DE gene counts per cell type.
    Positive bars (right): genes higher in group1 (CORT).
    Negative bars (left):  genes higher in group2 (OIL).

    FDR is computed per cell type via BH correction on (1 - proba_m1) / (1 - proba_m2).

    Parameters
    ----------
    input_dir  : directory containing de_<cell_type>.csv files
    output_path: output PNG path
    fdr_thresh : BH-corrected FDR threshold
    lfc_thresh : |lfc| threshold applied alongside FDR
    min_genes  : minimum total DE genes to include a cell type
    sort_by    : "total" | "up" | "down"
    figsize    : figure size
    group1     : label for the first group (positive/right bars)
    group2     : label for the second group (negative/left bars)
    """
    csv_files = sorted(Path(input_dir).glob("de_*.csv"))
    if not csv_files:
        print(f"No de_*.csv files found in {input_dir}")
        return

    summary = []
    for fpath in csv_files:
        try:
            df = pd.read_csv(fpath, index_col=0)
            required = {"proba_m1", "proba_m2", "lfc"}
            if not required.issubset(df.columns):
                print(f"  Skipping {fpath.name}: missing columns {required - set(df.columns)}")
                continue

            # BH FDR correction: treat (1 - proba) as local false-positive rate
            _, fdr_up,   _, _ = multipletests(1 - df["proba_m1"], method="fdr_bh")
            _, fdr_down, _, _ = multipletests(1 - df["proba_m2"], method="fdr_bh")

            up   = int(((fdr_up   < fdr_thresh) & (df["lfc"] >  lfc_thresh)).sum())
            down = int(((fdr_down < fdr_thresh) & (df["lfc"] < -lfc_thresh)).sum())

            if (up + down) >= min_genes:
                ct = fpath.stem.replace("de_", "", 1)
                summary.append({"cell_type": ct, "up": up, "down": down})

        except Exception as e:
            print(f"  Skipping {fpath.name}: {e}")

    if not summary:
        print("No cell types met the threshold.")
        return

    df_plot = pd.DataFrame(summary)
    df_plot["total"] = df_plot["up"] + df_plot["down"]

    sort_col = {"total": "total", "up": "up", "down": "down"}.get(sort_by, "total")
    df_plot = df_plot.sort_values(sort_col, ascending=True)

    fig, ax = plt.subplots(figsize=figsize)

    y = np.arange(len(df_plot))

    ax.barh(y,  df_plot["up"],           color="#C0392B", label=f"Higher in {group1}")
    ax.barh(y, -df_plot["down"].values,  color="#2980B9", label=f"Higher in {group2}")

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(df_plot["cell_type"], fontsize=6)
    ax.set_xlabel("Number of DE genes")
    ax.set_title(
        f"scANVI DE: {group1} vs {group2} (BH FDR < {fdr_thresh}, |lfc| > {lfc_thresh})",
        fontsize=10,
    )
    ax.legend(fontsize=8)

    # Mirror x-axis labels to show absolute counts on both sides
    xticks = ax.get_xticks()
    ax.set_xticklabels([str(int(abs(x))) for x in xticks])

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved: {output_path}")


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir",    "-i", default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output",       "-o", default=DEFAULT_OUTPUT)
    parser.add_argument("--fdr-thresh",   type=float, default=0.05)
    parser.add_argument("--lfc-thresh",   type=float, default=0.25)
    parser.add_argument("--min-genes",    type=int,   default=1)
    parser.add_argument("--sort-by",      choices=["total", "up", "down"], default="total")
    parser.add_argument("--group1",       default="CORT")
    parser.add_argument("--group2",       default="OIL")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    plot_scanvi_de_barplot(
        input_dir=args.input_dir,
        output_path=args.output,
        fdr_thresh=args.fdr_thresh,
        lfc_thresh=args.lfc_thresh,
        min_genes=args.min_genes,
        sort_by=args.sort_by,
        group1=args.group1,
        group2=args.group2,
    )
