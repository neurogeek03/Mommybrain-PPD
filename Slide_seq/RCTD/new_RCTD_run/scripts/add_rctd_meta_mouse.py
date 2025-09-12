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
import argparse

print('libraries loaded!')

# ========== ARGS ==========
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help="Path to input i.e merged rctd .csv files from the various subsets")
parser.add_argument("-o", "--output_dir", help="Path to output folder to store the anndata objects w/ metadata", default=None)
args = parser.parse_args()

# ========== PATHS ==========
project_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD'
data_folder = os.path.join(project_dir, 'genes_mouse_rename/slide_seq_full_adata_with_mouse_orthologs')
# rctd_out_folder = os.path.join(project_dir, 'new_RCTD_run/out_RCTD_all/test_ratified_ref_merged')
rctd_out_folder =args.input_dir
# output_dir = os.path.join(project_dir,'new_RCTD_run', 'ratified_ref_anndata_objects')
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

# ========== FIND SAMPLES ==========
# Identify all folders like B01_3TB, B08_3TB
# samples = [d.replace('_3TB', '') for d in os.listdir(data_folder) if d.endswith('_3TB')]
# samples = ['B37']

# samples = [fname.split('.')[0] for fname in os.listdir(data_folder)]
# ========== GET FILES ==========
rctd_csvs = sorted([f for f in os.listdir(rctd_out_folder) if f.endswith("_merged_RCTD.csv")])
samples = [f.split("_")[0] for f in rctd_csvs]
print(samples)

# samples = []
# # Loop through all files in the directory
# for filename in os.listdir(rctd_out_folder):
#     if os.path.isfile(os.path.join(rctd_out_folder, filename)):
#         parts = filename.split("_")
#         if len(parts) > 3:
#             sample = parts[3]  # 0-based index: 3 = 4th element
#             samples.append(sample)

# # Print the list of samples
# print(samples)

# ========== LOOP THROUGH SAMPLES ==========
for sample in samples:
    print(f"🔄 Processing sample: {sample}")

    h5ad_path = os.path.join(data_folder, f"{sample}.h5ad")
    rctd_out_path = os.path.join(rctd_out_folder, f"{sample}_merged_RCTD.csv")
    output_path = os.path.join(output_dir, f"{sample}_with_RCTD_mouse.h5ad")

    if not os.path.exists(h5ad_path):
        print(f"❌ h5ad not found for {sample}")
        continue

    if not os.path.exists(rctd_out_path):
        print(f"⚠️  RCTD file missing for {sample}, skipping RCTD annotation...")
        rctd_df = None
    else:
        rctd_df = pd.read_csv(rctd_out_path, index_col=0)

        print(rctd_df.head())
        #rctd_df = rctd_df.set_index("Unnamed: 0")
        print('rctd df columns:')
        print(rctd_df.columns)

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

    if rctd_df is not None:
        # Check match stats
        n_total = len(rctd_df)
        n_matched = rctd_df.index.isin(adata.obs.index).sum()
        print(f"🧬 RCTD: matched {n_matched} of {n_total} barcodes ({n_matched / n_total:.2%})")

        if n_matched == 0:
            print("⚠️  No barcodes matched! Skipping RCTD merge.")
        else:
            # Keep only the columns you want and add _mouse suffix
            cols_to_keep = ["spot_class", "first_type", "second_type", "min_score", "singlet_score"]
            rctd_df = rctd_df[cols_to_keep].rename(
                columns={col: "RCTD_" + col + "_rat" for col in cols_to_keep}
            )

            # Join with adata.obs
            adata.obs = adata.obs.join(rctd_df, how="left")

    # if rctd_df is not None:
    #     # Check match stats
    #     n_total = len(rctd_df)
    #     n_matched = rctd_df.index.isin(adata.obs.index).sum()
    #     print(f"🧬 RCTD: matched {n_matched} of {n_total} barcodes ({n_matched / n_total:.2%})")

    # if n_matched == 0:
    #     print("⚠️  No barcodes matched! Skipping RCTD merge.")
    # else:
    #     # Keep only the columns you want and add _mouse suffix
    #     cols_to_keep = ["spot_class", "first_type", "second_type", "min_score", "singlet_score"]
    #     rctd_df = rctd_df[cols_to_keep].rename(columns={col: "RCTD_" + col + "_rat" for col in cols_to_keep})

    #     # Join with adata.obs
    #     adata.obs = adata.obs.join(rctd_df, how="left")


    # ========== SAVE OUTPUT ==========
    adata.write(output_path)
    print(f"💾 Saved updated AnnData to: {output_path}\n")