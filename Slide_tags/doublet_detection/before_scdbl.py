"""
Title: Before scdbl finder - preparing anndata objects
Description:  Generates temporary anndata object files that contain raw count matrices
Author:   Maria Eleni Fafouti 
Date: 09-05-2025
"""
# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from rpy2.robjects import r, pandas2ri

# ========== PATHS ==========
ad_folder = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender'
out_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Doublet_detection/scanpy'

# sample_list = ["BC14"]
sample_list = [ "BC14", "BC28", "BC13", "BC15", "BC3", "BC9"]

for sample in sample_list:
    print(f'Sample {sample} is being processed')
    sample_folder = os.path.join(ad_folder, sample)
    h5ad_path = os.path.join(sample_folder, f"{sample}annotated_no_dim_red.h5ad")
    ad = sc.read_h5ad(h5ad_path)

    if ad.raw is not None:
        print("Switching .X to .raw.X (raw counts)")
        ad.X = ad.raw.X
        ad.var = ad.raw.var.copy()
    else:
        print(f"Warning: {sample} has no .raw slot!")

    print("Shape of .X:", ad.X.shape)
    print("Sparse matrix?" , issparse(ad.X))

    temp_path = os.path.join(out_dir, "tmp_dir", f"{sample}_tmp.h5ad")
    print(temp_path)
    ad.write_h5ad(temp_path)