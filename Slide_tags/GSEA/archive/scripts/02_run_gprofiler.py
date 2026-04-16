"""
Title: Run gProfiler enrichment per cell type
Description: Reads per-cell-type .rnk files (ordered mode) or edgeR results
             (ora mode) and runs gProfiler enrichment. Results saved as one
             CSV per cell type.
Author: Maria Eleni Fafouti
Date: 2026-04-02

Install:
  uv pip install gprofiler-official

Usage:
  # Ordered query (ranked list, all genes)
  python 02_run_gprofiler.py <rnk_dir> <output_dir> --mode ordered

  # ORA (unordered, significant DE genes only)
  python 02_run_gprofiler.py <rnk_dir> <output_dir> --mode ora --edger-dir <edger_results_dir>

Arguments:
  input_dir        directory containing per-cell-type .rnk files
  output_dir       directory where gProfiler CSVs will be saved
  --mode           ordered | ora (default: ordered)
  --edger-dir      directory with *_edgeR_results.tsv files (required for ora mode)
  --de-fdr         FDR threshold for selecting DE genes in ora mode (default: 0.05)
  --organism       gProfiler organism ID (default: hsapiens)
  --background-mode
    global         union of all genes detected across all cell types (default)
    per-celltype   each cell type uses only its own detected genes
    none           full genome (no background correction)
"""

import argparse
import pandas as pd
from pathlib import Path
from gprofiler import GProfiler

# ========== ARGS ==========
parser = argparse.ArgumentParser(description="Run gProfiler enrichment per cell type.")
parser.add_argument("input_dir", type=Path,
                    help="Directory containing per-cell-type .rnk files.")
parser.add_argument("output_dir", type=Path,
                    help="Directory where gProfiler result CSVs will be saved.")
parser.add_argument("--mode", choices=["ordered", "ora"], default="ordered",
                    help="ordered: ranked query; ora: unordered significant DE genes (default: ordered).")
parser.add_argument("--edger-dir", type=Path, default=None,
                    help="Directory with *_edgeR_results.tsv files (required for ora mode).")
parser.add_argument("--de-fdr", type=float, default=0.05,
                    help="edgeR FDR threshold for selecting DE genes in ora mode (default: 0.05).")
parser.add_argument("--organism", type=str, default="hsapiens",
                    help="gProfiler organism ID (default: hsapiens).")
parser.add_argument("--background-mode", choices=["global", "per-celltype", "none"],
                    default="global",
                    help="Background gene universe (default: global).")
parser.add_argument("--gprofiler-fdr", type=float, default=0.05,
                    help="gProfiler FDR threshold for enrichment significance (default: 0.05).")
args = parser.parse_args()

if args.mode == "ora" and args.edger_dir is None:
    parser.error("--edger-dir is required when --mode ora")

args.output_dir.mkdir(parents=True, exist_ok=True)

rnk_files = sorted(args.input_dir.glob("*.rnk"))
print(f"Found {len(rnk_files)} .rnk files in {args.input_dir}")
print(f"Mode: {args.mode} | Background: {args.background_mode}")

# ========== GLOBAL BACKGROUND ==========
if args.background_mode == "global":
    all_genes = set()
    for f in rnk_files:
        genes = pd.read_csv(f, sep="\t", header=None, index_col=0).index.tolist()
        all_genes.update(genes)
    global_background = list(all_genes)
    print(f"Global background: {len(global_background)} genes")

gp = GProfiler(return_dataframe=True)

for f in rnk_files:
    celltype = f.stem
    print(f"\nRunning gProfiler: {celltype}")

    all_tested_genes = pd.read_csv(f, sep="\t", header=None, index_col=0).index.tolist()

    # ---- Query genes ----
    if args.mode == "ora":
        edger_file = args.edger_dir / f"{celltype}_edgeR_results.tsv"
        if not edger_file.exists():
            print(f"  Skipping — edgeR file not found: {edger_file.name}")
            continue
        de_df = pd.read_csv(edger_file, sep="\t", index_col=0)
        query_genes = de_df[de_df["FDR"] < args.de_fdr].index.tolist()
        print(f"  DE genes (FDR < {args.de_fdr}): {len(query_genes)}")
        if len(query_genes) == 0:
            print(f"  Skipping — no DE genes at threshold")
            continue
    else:
        query_genes = all_tested_genes

    # ---- Background ----
    if args.background_mode == "global":
        background = global_background
    elif args.background_mode == "per-celltype":
        background = all_tested_genes
    else:
        background = None

    results = gp.profile(
        organism=args.organism,
        query=query_genes,
        background=background,
        ordered=(args.mode == "ordered"),
        significance_threshold_method="fdr",
        user_threshold=args.gprofiler_fdr,
        no_evidences=False,
    )

    if results.empty:
        print(f"  No significant terms for {celltype}")
        continue

    out_path = args.output_dir / f"{celltype}_gprofiler.csv"
    results.to_csv(out_path, index=False)
    print(f"  {len(results)} terms. Saved: {out_path.name}")

print(f"\nDone. Results saved to: {args.output_dir}")
