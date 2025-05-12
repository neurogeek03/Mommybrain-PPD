"""
Title: Before mast - preparing anndata object
Description:  Generates a normalized, log transformed (not scaled) anndata object with a complete var slot (gene symbols, etc)
Author:   Maria Eleni Fafouti 
Date: 11-05-2025
"""
# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import scanpy as sc
import os 
import scipy.sparse
import numpy as np
import pandas as pd

# ========== PATHS ==========
adata_path_target = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Doublet_detection/scanpy/doublets_harmony.h5ad"
out_dir = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/DE/MAST"

# ========== READING & PRE-PROCESSING DATA ==========
adata = sc.read_h5ad(adata_path_target)

# Copying raw to X 
adata.X = adata.raw.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy() # saving as new layer

# print(adata.shape)  # should be (207903, 25629)
# print(len(adata.obs_names))  # should be 207903
# print(len(adata.var_names)) 

# Filter out doublets
adata_filtered = adata[adata.obs['scDblFinder.class'] == "singlet"].copy()

print(adata_filtered.var.head())


# EXTRACTING VAR 
adata_path = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender/BC9/converted_ann_data_BC9.h5ad"

adata_older = sc.read_h5ad(adata_path)
print(adata_older.var.head())

adata_filtered.var = adata_older.var.copy()
print(adata_filtered.var.head())

out_path = os.path.join(out_dir, "var_noscale_harmony_doubletfiltered.h5ad")
adata_filtered.write(out_path)