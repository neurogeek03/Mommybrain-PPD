"""
make_liana_dea_csv.py
---------------------
Converts per-cell-type EdgeR LRT results into the LIANA-compatible DEA CSV
format (same schema as out/first_try/dea_result.csv produced by PyDESeq2).

Column mapping:
  EdgeR            -> LIANA/DESeq2 schema
  -------            ----------------------
  (row index)      -> gene (index)
  (filename stem)  -> subclass_name  [spaces restored from underscores]
  logFC            -> log2FoldChange and stat  [logFC used as signed stat for LIANA directionality]
  logCPM           -> baseMean        [proxy; different scale — not used by df_to_lr]
  LR               -> (not used)      [LR is unsigned; replaced by logFC as stat]
  PValue           -> pvalue          [used by df_to_lr]
  FDR              -> padj            [used by df_to_lr]
  (absent)         -> lfcSE = NaN     [SE of log2FC; not used by df_to_lr]

LIANA df_to_lr uses stat_keys=['stat', 'pvalue', 'padj'] (CCC_liana.py:208).
stat=logFC gives interaction_stat a direction (positive=up in CORT, negative=down in CORT),
matching the signed Wald stat used by PyDESeq2.

Usage:
  python make_liana_dea_csv.py <dea_dir>

  dea_dir  path to directory containing *_edgeR_results.tsv files
           (e.g. ../out/edger_lrt or ../out/new_march_26/edger_lrt)

Output is written to <dea_dir>/../../edgeR_dge_input_liana.csv
(i.e. alongside the 'out/' directory, relative to dea_dir).

Requires: conda activate liana
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser(description="Convert EdgeR LRT results to LIANA-compatible DEA CSV.")
parser.add_argument("dea_dir", type=Path, help="Directory containing *_edgeR_results.tsv files")
args = parser.parse_args()

dea_dir  = args.dea_dir.resolve()
out_path = dea_dir.parent / "edgeR_dge_input_liana.csv"

tsv_files = sorted(dea_dir.glob("*_edgeR_results.tsv"))
print(f"Found {len(tsv_files)} cell-type files in {dea_dir}")

frames = []
for f in tsv_files:
    # Recover original subclass_name: filename used spaces→underscores
    ct = f.name.replace("_edgeR_results.tsv", "").replace("_", " ")

    df = pd.read_csv(f, sep="\t", index_col=0)
    df.index.name = None  # gene name is the index

    out_df = pd.DataFrame({
        "subclass_name":   ct,
        "baseMean":        df["logCPM"],          # log-CPM; proxy for baseMean
        "log2FoldChange":  df["logFC"],
        "lfcSE":           np.nan,                # no equivalent in EdgeR LRT
        "stat":            df["logFC"],             # signed log2FC for LIANA directionality
        "pvalue":          df["PValue"],
        "padj":            df["FDR"],
    }, index=df.index)

    frames.append(out_df)

result = pd.concat(frames)
result.index.name = "index"

result.to_csv(out_path, index=True)
print(f"Saved {len(result):,} rows ({result['subclass_name'].nunique()} cell types) to:\n  {out_path}")
print(f"\nColumn summary:\n{result.dtypes}")
print(f"\nFirst rows:\n{result.head(3)}")
