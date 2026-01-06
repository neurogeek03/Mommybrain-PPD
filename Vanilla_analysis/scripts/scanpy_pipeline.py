import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
data_path = Path('/scratch/mfafouti/Mommybrain/TIA_toolbox/out/FINAL/objects')
figures_dir = output_base / 'figures'
figures_dir.mkdir(exist_ok=True, parents=True)

# =================== PARAMS ===================
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor="white")
sc.settings.figdir = figures_dir
res=1
min_genes=50

sample = 'B01'

# =================== INPUT ===================
sample_files = list(data_path.glob(f"{sample}*"))
adata_path = sample_files[0]
adata = sc.read_h5ad(adata_path)

# sanity
adata 

# # sanity: plot spatial to check
# fig, ax = plt.subplots(figsize=(5, 5))
# x = adata_tissue.obsm["X_spatial"][:, 0]
# y = adata_tissue.obsm["X_spatial"][:, 1]
# ax.scatter(x, y, s=1)
# ax.invert_yaxis()
# ax.set_aspect("equal")
# png1_path = output_base / f"{sample}_test.png"
# fig.savefig(png1_path, dpi=300)

# =================== Add metadata ===================
# clean gene names and IDs
split = adata.var_names.to_series().str.split("-", n=1, expand=True)
adata.var_names = split[1] # keep second part (ensembl ID)
adata.var["name"] = split[0].values # keep name only (if absent, theres a ".")
adata
adata.var.head()

# label mitochondrial genes 
adata.var["mt_gene"] = adata.var["name"].str.upper().str.startswith("Mt-")
adata.var.head()

# # sanity
# mt_mask = adata.var["name"].str.startswith("MT-")
# print(adata.var[mt_mask].head())

# intermediate save 
adata_path = output_base / 'objects' / f'{sample}_all_metadata.h5ad'
adata.write_h5ad(adata_path)

# =================== QC ===================
# SUBSET to on-tissue
adata_tissue = adata[adata.obs["in_tissue"]].copy()

adata_tissue.obs["nCount_RNA"] = adata_tissue.X.sum(axis=1).A1
adata_tissue.obs["nFeature_RNA"] = (adata_tissue.X > 0).sum(axis=1).A1

# Mitochondrial genes (human convention: MT-)
sc.pp.calculate_qc_metrics(
    adata_tissue,
    qc_vars=["mt_gene"],
    percent_top=None,
    log1p=False,
    inplace=True
)

sc.pl.violin(
    adata_tissue,
    keys=["total_counts", "n_genes_by_counts", "pct_counts_mt_gene"],
    jitter=0.4,
    show=False,
    multi_panel=True,
    save="_qc_violin.png"
)

# Filter cells and genes
sc.pp.filter_cells(adata_tissue, min_genes=min_genes)
sc.pp.filter_genes(adata_tissue, min_cells=3)
adata_tissue = adata_tissue[adata_tissue.obs.pct_counts_mt_gene < 10].copy()

# Normalization and log transform
sc.pp.normalize_total(adata_tissue, target_sum=1e4)
sc.pp.log1p(adata_tissue)

# Highly variable genes
sc.pp.highly_variable_genes(
    adata_tissue,
    flavor="seurat_v3",
    n_top_genes=2000
)
adata_tissue = adata_tissue[:, adata_tissue.var.highly_variable].copy()

# Scaling and PCA
sc.pp.scale(adata_tissue, max_value=10)
sc.tl.pca(adata_tissue, svd_solver="arpack")

# Neighborhood graph
sc.pp.neighbors(adata_tissue, n_neighbors=10, n_pcs=40)

# Embedding
sc.tl.umap(adata_tissue)

# Clustering
sc.tl.leiden(adata_tissue, resolution=res)

# Marker gene detection
sc.tl.rank_genes_groups(
    adata_tissue,
    groupby="leiden",
    method="wilcoxon"
)

# Plotting
sc.pl.umap(adata_tissue, color=["leiden"], save="res_1_mingenes_50_leiden.png")
sc.pl.rank_genes_groups(adata_tissue, n_genes=20, sharey=False, save="_rank_genes_groups.png")

# Save processed object
ad_proc_path = output_base / 'objects' / f"{sample}_res_{res}_mingenes_{min_genes}_adata_processed.h5ad"
adata_tissue.write(ad_proc_path)
