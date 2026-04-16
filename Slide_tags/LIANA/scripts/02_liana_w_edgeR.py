import argparse
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import glob

import liana as li

parser = argparse.ArgumentParser(description='Run LIANA df_to_lr with EdgeR DEA results')
parser.add_argument('--dea_path',      type=str, required=True,  help='Path to EdgeR DEA CSV')
parser.add_argument('--celltype_key',  type=str, default='subclass_name')
parser.add_argument('--condition_key', type=str, default='group')
parser.add_argument('--resource_name', type=str, default='mouseconsensus')
parser.add_argument('--expr_prop',     type=float, default=0.1)
parser.add_argument('--run_dir',       type=str,   required=True,
                    help='Path to the run output directory (e.g. out/runs/my_run)')
args = parser.parse_args()

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"Current working directory: {project_path}")
out_dir  = Path(args.run_dir) / 'liana_edgeR'
out_dir.mkdir(exist_ok=True, parents=True)
figs_dir = out_dir / 'figures'
figs_dir.mkdir(exist_ok=True, parents=True)
sc.settings.figdir = figs_dir

# =================== INPUT ===================
# Preprocessed AnnData from 00_preprocess_adata.py
preprocess_dir = Path(args.run_dir) / 'preprocessing'
h5ad_files = glob.glob(str(preprocess_dir / 'object_LIANA_*.h5ad'))
if len(h5ad_files) != 1:
    raise ValueError(f"Expected exactly 1 preprocessed .h5ad in {preprocess_dir}, found: {h5ad_files}")
adata_path = Path(h5ad_files[0])
print(f'Using {adata_path} ...')
adata = sc.read_h5ad(adata_path)

# EdgeR DEA results
dea_df = pd.read_csv(args.dea_path, index_col=0)
print(f'Using DEA from {args.dea_path} ({len(dea_df)} rows)')

# =================== LIANA df_to_lr ===================
lr_res = li.multi.df_to_lr(
    adata,
    dea_df=dea_df,
    resource_name=args.resource_name,  # NOTE: uses MOUSE gene symbols
    expr_prop=args.expr_prop,           # expression proportion filter (pooled across conditions — see CLAUDE.md)
    groupby=args.celltype_key,
    stat_keys=['stat', 'pvalue', 'padj'],
    use_raw=False,
    complex_col='stat',                 # Wald stat used to resolve complexes
    verbose=True,
    return_all_lrs=False,
)

lr_res = lr_res.sort_values("interaction_stat", ascending=False, key=abs)

lr_res_path = out_dir / 'lr_res.csv'
lr_res.to_csv(lr_res_path)
print(f"Saved lr_res to {lr_res_path}")
