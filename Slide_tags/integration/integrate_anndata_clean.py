"""
Title: Scanpy Integration  
Description:  Integrating the anndata objects from multiple samples 
Author:   Maria Eleni Fafouti 
Date: 07-05-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import anndata
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import harmonypy as hm
import bbknn
import numpy as np

# ========== DEFINING ESSENTIAL PATHS ==========
ad_folder = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender'
out_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Integration/scanpy'
csv_path = os.path.join(ad_folder, "hex_codes_subclass.csv")

df = pd.read_csv(csv_path, header=None, names=['cell_type', 'color'])
color_dict = dict(zip(df['cell_type'], df['color']))

sample_list = ['BC28', 'BC3', 'BC9', 'BC15', 'BC14', 'BC13']
# sample_list = ['BC28', 'BC3']

# ========== CONCANTENATING ADATA FILES ==========
adata_list = []

for sample in sample_list: 
    print(f"Processing {sample}...")

    # Loading anndata objects 
    sample_folder = os.path.join(ad_folder, sample)
    h5ad_path = os.path.join(sample_folder, f"{sample}annotated_no_dim_red.h5ad")
    ad = sc.read_h5ad(h5ad_path)

    print("X shape:")
    print(ad.X.shape)  # we will keep the raw matric in each anndata obj
    print("Raw matrix:")
    print(ad.raw.X.shape)
    print("Observation metadata (obs):")
    print(ad.obs.columns.tolist())

    # If .raw exists, restore .X from .raw.X and update .var accordingly
    if ad.raw is not None:
        print("Switching .X to .raw.X (raw counts)")
        ad.X = ad.raw.X
        ad.var = ad.raw.var.copy()
    else:
        print(f"Warning: {sample} has no .raw slot!")

    # Save to list 
    adata_list.append(ad)

# Concatenating anndata objects 
adata_all = anndata.concat(adata_list)

# ========== FINDING HIGHLY VARIABLE GENES (HVG) ==========
# Saving the merged data as the raw matrix first 
adata_all.raw = adata_all.copy() 

# Batch-aware HVG selection
# flavor='seurat_v3' => the only one which is multi-sample and batch aware
sc.pp.highly_variable_genes(adata_all.raw.to_adata(), batch_key='sample', flavor='seurat_v3', n_top_genes=2000)

# ========== NORMALIZATION AND LOG TRANSFORM ==========
# Normalize total counts per cell (so every cell has the same library size)
sc.pp.normalize_total(adata_all, target_sum=1e4)

# Log-transform the data (this helps reduce variance due to outliers)
sc.pp.log1p(adata_all)

# ========== SCALING AND PCA ==========
sc.pp.scale(adata_all, max_value=10)
sc.tl.pca(adata_all, svd_solver="arpack")

# # ========== RUNNING HARMONY ==========
# sce.pp.harmony_integrate(adata_all, key="sample", max_iter_harmony=20)  # Creates `X_pca_harmony` in adata.obsm

# # ========== USING HARMONY-ADJUSTED PCs DOWNSTREAM: NEIGHBORS & CLUSTERING ==========
# sc.pp.neighbors(adata_all, use_rep="X_pca_harmony")
# sc.tl.umap(adata_all)
# sc.tl.leiden(adata_all)


# ========== RUNNING BBKNN ==========
sc.external.pp.bbknn(adata_all, batch_key="sample", computation="full") # computation ="approx" for faster neighborhood approximations
sc.tl.umap(adata_all)
sc.tl.leiden(adata_all) 

new_file_dir = os.path.join(out_dir, "integrated_bbknn.h5ad")
adata_all.write(new_file_dir)

sc.pl.umap(adata_all, color=['MapMyCells_cell_type'], palette=color_dict)
plt.suptitle("BBKNN UMAPs", fontsize=16)
fig_path = os.path.join(out_dir, "umap_BBKNN_all_colors.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', pad_inches=0.1)

# # ========== PLOTTING ==========
# sc.pl.umap(adata_all, color='MapMyCells_cell_type', palette=color_dict)
# plt.title(f"Harmony - color: MapMyCells", fontsize=14)
# fig_path = os.path.join(out_dir, "umap_Harmony_mapmycells.png")
# plt.savefig(fig_path, dpi=300, bbox_inches='tight', pad_inches=0.1)

# sc.pl.umap(adata_all, color='sample', palette=color_dict)
# plt.title(f"Harmony - color: sample ID", fontsize=14)
# fig_path = os.path.join(out_dir, "umap_Harmon_sampleID.png")
# plt.savefig(fig_path, dpi=300, bbox_inches='tight', pad_inches=0.1)

# sc.pl.umap(adata_all, color='treatment', palette=color_dict)
# plt.title(f"Harmony - color: treatment", fontsize=14)
# fig_path = os.path.join(out_dir, "umap_Harmony_treatment.png")
# plt.savefig(fig_path, dpi=300, bbox_inches='tight', pad_inches=0.1)