"""
Title: Preview slide-seq anndata objects & add metadata
Description:  Add RCTD info (i.e cell types)
Author:   Maria Eleni Fafouti 
Date: 30-06-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import os
import anndata as ad
import pandas as pd
import numpy as np

# ========== PATHS ==========
data_folder = '/project/s/shreejoy/MommyBrain/CurioSeeker_out_rat'
project_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Spatial'
rctd_out_folder = os.path.join(data_folder, 'RCTD_out_merged')
output_dir = os.path.join(project_dir, 'annadata_objects')
os.makedirs(output_dir, exist_ok=True)

# ========== FIND SAMPLES ==========
# Identify all folders like B01_3TB, B08_3TB
# samples = [d.replace('_3TB', '') for d in os.listdir(data_folder) if d.endswith('_3TB')]
samples = ['B41']

# ========== LOOP THROUGH SAMPLES ==========
for sample in samples:
    print(f"🔄 Processing sample: {sample}")

    h5ad_path = os.path.join(data_folder, f"{sample}_3TB", sample, f"{sample}_anndata.h5ad")
    rctd_out_path = os.path.join(rctd_out_folder, f"{sample}_merged_RCTD.csv")
    output_path = os.path.join(output_dir, f"{sample}_with_RCTD.h5ad")

    if not os.path.exists(h5ad_path):
        print(f"❌ h5ad not found for {sample}")
        continue

    if not os.path.exists(rctd_out_path):
        print(f"⚠️  RCTD file missing for {sample}, skipping RCTD annotation...")
        rctd_df = None
    else:
        rctd_df = pd.read_csv(rctd_out_path, index_col=0)

    # ========== READ H5AD ==========
    adata = sc.read_h5ad(h5ad_path)
    print(f"✅ Loaded {h5ad_path} with shape {adata.shape}")

    # ========== CLEAN & SPLIT VAR NAMES ==========
    split_info = adata.var_names.str.extract(r'^(?:(.*)-)?(ENSRNOG\d+)$')
    split_info[0] = split_info[0].fillna(split_info[1])
    split_info.columns = ['gene_symbol', 'gene_id']
    adata.var['gene_symbol'] = split_info['gene_symbol'].values
    adata.var['gene_id'] = split_info['gene_id'].values
    adata.var_names = adata.var['gene_id']
    adata.var_names_make_unique()

    # ========== ADD RCTD METADATA ==========
    if rctd_df is not None:
        # Check match stats
        n_total = len(rctd_df)
        n_matched = rctd_df.index.isin(adata.obs.index).sum()
        print(f"🧬 RCTD: matched {n_matched} of {n_total} barcodes ({n_matched / n_total:.2%})")

        if n_matched == 0:
            print("⚠️  No barcodes matched! Skipping RCTD merge.")
        else:
            # Clean and merge
            rctd_df = rctd_df[["spot_class", "first_type", "second_type", "min_score", "singlet_score"]]
            rctd_df.columns = ["RCTD_" + col for col in rctd_df.columns]
            adata.obs = adata.obs.join(rctd_df, how="left")

    # ========== SAVE OUTPUT ==========
    adata.write(output_path)
    print(f"💾 Saved updated AnnData to: {output_path}\n")
