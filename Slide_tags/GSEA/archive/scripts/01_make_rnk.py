"""
Title: Make pre-ranked gene lists for GSEA
Description: Reads per-cell-type EdgeR LRT results and writes one .rnk file
             per cell type, ranked by sign(logFC) * -log10(PValue).
Author: Maria Eleni Fafouti
Date: 2026-04-02

Usage:
  python 06_make_rnk.py <input_dir> <output_dir>

  input_dir   directory containing *_edgeR_results.tsv files
              (e.g. runs/merge_imn/edger_lrt/OIL_vs_CORT)
  output_dir  directory where per-cell-type .rnk files will be saved
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# ========== ARGS ==========
parser = argparse.ArgumentParser(description="Make pre-ranked gene lists for GSEA from EdgeR results.")
parser.add_argument("input_dir", type=Path,
                    help="Directory containing *_edgeR_results.tsv files.")
parser.add_argument("output_dir", type=Path,
                    help="Directory where per-cell-type .rnk files will be saved.")
args = parser.parse_args()

args.output_dir.mkdir(parents=True, exist_ok=True)

tsv_files = sorted(args.input_dir.glob("*_edgeR_results.tsv"))
print(f"Found {len(tsv_files)} cell-type files in {args.input_dir}")


for f in tsv_files:
    celltype = f.name.replace("_edgeR_results.tsv", "")

    df = pd.read_csv(f, sep="\t", index_col=0)

    # Ranking metric: sign(logFC) * -log10(PValue)
    df["rank_metric"] = np.sign(df["logFC"]) * -np.log10(df["FDR"])
    df = df[["rank_metric"]].sort_values("rank_metric", ascending=False)
    df.index.name = None

    out_path = args.output_dir / f"{celltype}.rnk"
    df.to_csv(out_path, sep="\t", header=False)
    print(f"  Written: {out_path.name}  ({len(df)} genes)")

print(f"\nDone. .rnk files saved to: {args.output_dir}")
