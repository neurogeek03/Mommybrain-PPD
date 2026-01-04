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

# =================== PARAMS ===================
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor="white")

sample = 'B01'

# =================== INPUT ===================
sample_files = list(data_path.glob(f"{sample}*"))
adata_path = sample_files[0]

adata = sc.read_h5ad(adata_path)
adata_tissue = adata[adata.obs["in_tissue"]].copy()
fig, ax = plt.subplots(figsize=(5, 5))

x = adata_tissue.obsm["X_spatial"][:, 0]
y = adata_tissue.obsm["X_spatial"][:, 1]

ax.scatter(x, y, s=1)
ax.invert_yaxis()
ax.set_aspect("equal")

png1_path = output_base / f"{sample}_test.png"
fig.savefig(png1_path, dpi=300)
# # =================== QC ===================

# # SUBSET to on-tissue


# adata.obs["n_counts"] = adata.X.sum(axis=1).A1
# adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1

# # Mitochondrial genes (human convention: MT-)
# adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
# sc.pp.calculate_qc_metrics(
#     adata,
#     qc_vars=["mt"],
#     percent_top=None,
#     log1p=False,
#     inplace=True
# )

# # Filter cells and genes
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# adata = adata[adata.obs.pct_counts_mt < 10].copy()

# # ----------------------------
# # Normalization and log transform
# # ----------------------------
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# # ----------------------------
# # Highly variable genes
# # ----------------------------
# sc.pp.highly_variable_genes(
#     adata,
#     flavor="seurat_v3",
#     n_top_genes=2000
# )
# adata = adata[:, adata.var.highly_variable].copy()

# # ----------------------------
# # Scaling and PCA
# # ----------------------------
# sc.pp.scale(adata, max_value=10)
# sc.tl.pca(adata, svd_solver="arpack")

# # ----------------------------
# # Neighborhood graph
# # ----------------------------
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# # ----------------------------
# # Embedding
# # ----------------------------
# sc.tl.umap(adata)

# # ----------------------------
# # Clustering
# # ----------------------------
# sc.tl.leiden(adata, resolution=0.5)

# # ----------------------------
# # Marker gene detection
# # ----------------------------
# sc.tl.rank_genes_groups(
#     adata,
#     groupby="leiden",
#     method="wilcoxon"
# )

# # ----------------------------
# # Plotting
# # ----------------------------
# sc.pl.umap(adata, color=["leiden"], save="_leiden.png")
# sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# # ----------------------------
# # Save processed object
# # ----------------------------
# adata.write("adata_processed.h5ad")