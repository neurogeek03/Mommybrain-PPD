import argparse
from datetime import datetime
from pathlib import Path
import scanpy as sc

parser = argparse.ArgumentParser(description='Preprocess AnnData for LIANA CCC pipeline')
parser.add_argument('--cort_samples',    nargs='+', required=True)
parser.add_argument('--oil_samples',     nargs='+', required=True)
parser.add_argument('--sample_key',      type=str,  default='sample')
parser.add_argument('--celltype_key',    type=str,  default='subclass_name')
parser.add_argument('--condition_key',   type=str,  default='group')
parser.add_argument('--n_top_genes',     type=int,  default=2000)
parser.add_argument('--merge_celltypes', nargs='+', default=None,
                    help='Cell type labels to merge')
parser.add_argument('--merge_label',     type=str,  default=None,
                    help='Label for merged cell types (default: joined with "+")')
parser.add_argument('--adata_path',      type=str,  required=True,
                    help='Path to input .h5ad file')
parser.add_argument('--run_dir',         type=str,  required=True,
                    help='Path to the run output directory (e.g. out/runs/my_run)')
args = parser.parse_args()

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"Current working directory: {project_path}")
data_dir = project_path / 'data'
out_dir  = Path(args.run_dir) / 'preprocessing'
out_dir.mkdir(exist_ok=True, parents=True)
figs_dir = out_dir / 'figures'
figs_dir.mkdir(exist_ok=True, parents=True)
sc.settings.figdir = figs_dir

# =================== INPUT ===================
adata_path = Path(args.adata_path)
adata = sc.read_h5ad(adata_path)
print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# =================== GROUP LABELS ===================
adata.obs[args.condition_key] = None
adata.obs.loc[adata.obs[args.sample_key].isin(args.cort_samples), args.condition_key] = 'CORT'
adata.obs.loc[adata.obs[args.sample_key].isin(args.oil_samples),  args.condition_key] = 'OIL'

unlabelled = adata.obs[args.condition_key].isna().sum()
if unlabelled > 0:
    print(f"WARNING: {unlabelled} cells could not be assigned a group label.")

# =================== MERGE CELL TYPES (optional) ===================
if args.merge_celltypes:
    merge_label = args.merge_label or "+".join(args.merge_celltypes)
    mask = adata.obs[args.celltype_key].isin(args.merge_celltypes)
    print(f"Merging {mask.sum()} cells from {args.merge_celltypes} into '{merge_label}'.")
    adata.obs[args.celltype_key] = adata.obs[args.celltype_key].astype(str)
    adata.obs.loc[mask, args.celltype_key] = merge_label

# =================== SET X TO RAW COUNTS ===================
adata.X = adata.layers["counts"]
adata.layers["counts"] = adata.X.copy()

# =================== NORMALIZATION ===================
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# =================== HVG ===================
sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, batch_key=args.sample_key)
sc.pl.highly_variable_genes(adata, save='_hvg.png')

# =================== PCA / NEIGHBORS / UMAP ===================
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color=[args.condition_key, args.sample_key, 'class_name', args.celltype_key],
    frameon=False,
    ncols=2,
    save='_preprocess.png'
)

# =================== OUTPUT ===================
out_path = out_dir / f'object_LIANA_{adata.n_obs}_pseudobulk_gene_names_var_condition.h5ad'
adata.write(out_path)
print(f"Saved preprocessed object to {out_path}")
