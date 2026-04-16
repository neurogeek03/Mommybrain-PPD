"""
Title: Plot gProfiler enrichment results
Description: Reads per-cell-type gProfiler CSVs, filters to GO:BP, and produces
             a clustered heatmap and dot plot for three groups:
               - NN   (non-neuronal, names ending in _NN)
               - Gaba (GABAergic, names ending in _Gaba)
               - Rest (Glut + IMN)
Author: Maria Eleni Fafouti
Date: 2026-04-02

Usage:
  python 08_plot_gprofiler.py <input_dir> <output_dir> [--top-n 30] [--min-celltypes 2]
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ========== ARGS ==========
parser = argparse.ArgumentParser(description="Plot gProfiler GO:BP results per cell type group.")
parser.add_argument("input_dir", type=Path,
                    help="Directory containing per-cell-type *_gprofiler.csv files.")
parser.add_argument("output_dir", type=Path,
                    help="Directory where plots will be saved.")
parser.add_argument("--top-n", type=int, default=30,
                    help="Top N terms to show per group (by min p-value). Default: 30.")
parser.add_argument("--min-celltypes", type=int, default=2,
                    help="Only show terms significant in at least this many cell types. Default: 2.")
parser.add_argument("--celltypes-csv", type=Path, default=None,
                    help="Optional CSV with a 'celltype' column. Only those cell types will be plotted.")
parser.add_argument("--source", type=str, default="GO:BP",
                    help="gProfiler source to filter (e.g. GO:BP, GO:MF, GO:CC, KEGG, REAC). Default: GO:BP.")
args = parser.parse_args()

args.output_dir.mkdir(parents=True, exist_ok=True)

# ========== LOAD ==========
csv_files = sorted(args.input_dir.glob("*_gprofiler.csv"))
print(f"Found {len(csv_files)} cell-type files.")

frames = []
for f in csv_files:
    celltype = f.stem.replace("_gprofiler", "")
    df = pd.read_csv(f)
    df = df[df["source"] == args.source][["name", "p_value", "intersection_size"]].copy()
    df["celltype"] = celltype
    frames.append(df)

if not frames:
    print(f"No {args.source} results found.")
    exit(0)

all_df = pd.concat(frames, ignore_index=True)

if args.celltypes_csv is not None:
    keep = pd.read_csv(args.celltypes_csv)["celltype"].tolist()
    found = all_df["celltype"].unique().tolist()
    matched = [ct for ct in keep if ct in found]
    unmatched = [ct for ct in keep if ct not in found]
    print(f"  CSV contains {len(keep)} cell types, {len(matched)} matched, {len(unmatched)} unmatched")
    if unmatched:
        print(f"  Unmatched: {unmatched[:5]}")
    if found:
        print(f"  Example gprofiler names: {found[:3]}")
    all_df = all_df[all_df["celltype"].isin(keep)]
    print(f"Filtered to {all_df['celltype'].nunique()} cell types from {args.celltypes_csv.name}")

# ========== ASSIGN GROUPS ==========
def assign_group(ct):
    if ct.endswith("_NN"):
        return "NN"
    elif ct.endswith("_Gaba"):
        return "Gaba"
    else:
        return "Rest"

all_df["group"] = all_df["celltype"].apply(assign_group)

# ========== PLOT FUNCTIONS ==========
def make_heatmap(log_pivot, title, out_path):
    log_pivot = log_pivot[sorted(log_pivot.columns)]
    row_h = max(8, len(log_pivot) * 0.3)
    col_w = max(10, len(log_pivot.columns) * 0.7)
    g = sns.clustermap(
        log_pivot,
        figsize=(col_w, row_h),
        cmap="YlOrRd",
        col_cluster=False,
        linewidths=0.3,
        linecolor="lightgrey",
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "-log10(p-value)"},
        dendrogram_ratio=(0.1, 0.15),
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis="y", labelsize=8)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=8)
    plt.suptitle(title, y=1.01, fontsize=12)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path.name}")


def make_dotplot(dot_df, title, out_path):
    celltypes = sorted(dot_df["celltype"].unique())
    dot_df = dot_df.copy()
    dot_df["celltype"] = pd.Categorical(dot_df["celltype"], categories=celltypes, ordered=True)
    dot_df = dot_df.sort_values("celltype")
    terms = dot_df["name"].unique()
    col_w = max(6, len(celltypes) * 0.7)
    row_h = max(8, len(terms) * 0.3)
    fig, ax = plt.subplots(figsize=(col_w, row_h))
    scatter = ax.scatter(
        dot_df["celltype"],
        dot_df["name"],
        c=dot_df["neg_log10_p"],
        s=dot_df["intersection"] * 0.5,
        cmap="YlOrRd",
        alpha=0.85,
        edgecolors="grey",
        linewidths=0.3,
    )
    plt.colorbar(scatter, ax=ax, label="-log10(p-value)")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_title(title, fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path.name}")


# ========== PLOT PER GROUP ==========
for group in ["NN", "Gaba", "Rest"]:
    group_df = all_df[all_df["group"] == group]
    if group_df.empty:
        print(f"\nNo data for group: {group}, skipping.")
        continue

    print(f"\nProcessing group: {group} ({group_df['celltype'].nunique()} cell types)")

    pivot = group_df.pivot_table(
        index="name", columns="celltype", values="p_value", aggfunc="min"
    )
    log_pivot = -np.log10(pivot).fillna(0)

    # Filter: significant in >= min_celltypes
    sig_counts = (pivot < 0.05).sum(axis=1)
    log_pivot = log_pivot[sig_counts >= args.min_celltypes]

    if log_pivot.empty:
        print(f"  No terms pass min_celltypes={args.min_celltypes} filter, skipping.")
        continue

    # Top N terms
    top_terms = log_pivot.max(axis=1).nlargest(args.top_n).index
    log_pivot = log_pivot.loc[top_terms]

    print(f"  {len(log_pivot)} terms x {len(log_pivot.columns)} cell types")

    source_label = args.source.replace(":", "_")
    title = f"{args.source} Enrichment — {group}"
    make_heatmap(log_pivot, title, args.output_dir / f"gprofiler_{source_label}_{group}_heatmap.png")

    # Dot plot data
    pivot_int = group_df.pivot_table(
        index="name", columns="celltype", values="intersection_size", aggfunc="max"
    )
    pivot_int = pivot_int.reindex(index=top_terms).fillna(0)
    melt_log = log_pivot.reset_index().melt(id_vars="name", var_name="celltype", value_name="neg_log10_p")
    melt_int = pivot_int.reset_index().melt(id_vars="name", var_name="celltype", value_name="intersection")
    dot_df = melt_log.merge(melt_int, on=["name", "celltype"])
    dot_df = dot_df[dot_df["neg_log10_p"] > 0]

    make_dotplot(dot_df, title, args.output_dir / f"gprofiler_{source_label}_{group}_dotplot.png")

print(f"\nDone. Plots saved to: {args.output_dir}")
