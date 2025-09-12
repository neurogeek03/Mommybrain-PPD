#!/usr/bin/env python

import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse

# ========== ARGS ==========
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help="Path to input annotated anndata objects w/ RCTD output")
parser.add_argument("-o", "--output_dir", help="Path to output folder for the figures", default=None)
args = parser.parse_args()

# ========== PATHS ==========
adata_dir = args.input_dir
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

# ========== GET FILES ==========
h5ad_files = sorted([f for f in os.listdir(adata_dir) if f.endswith("_with_RCTD_mouse.h5ad")])
spot_class_column = "RCTD_spot_class_rat"
type_column = "RCTD_first_type_rat"

print(f"Found {len(h5ad_files)} h5ad files")

# ========== COLOR PALETTE - ABC Atlas ==========
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)
# Create mapping from label number (prefix before _) to hex color
def get_num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except:
        return -1
color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)
# Sort CSV by numeric prefix
color_df = color_df.sort_values("num_prefix")
# Build dictionary mapping label to hex
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

# ========== LOOP OVER FILES ==========
for h5ad_file in h5ad_files:
    sample = h5ad_file.split("_")[0]
    print(f"📂 Reading: {sample}")
    adata_path = os.path.join(adata_dir, h5ad_file)
    adata = sc.read_h5ad(adata_path)
    # Keep only singlets
    if spot_class_column in adata.obs:
        adata = adata[adata.obs[spot_class_column] == "singlet"].copy()
    else:
        print(f"⚠️ {spot_class_column} not found in {sample}, skipping...")
        continue

    print("➡️ Starting computations...")
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    print(f"✅ Done with UMAP for {sample}")

    # Saving object
    adata_path = os.path.join(output_dir, "objects")
    os.makedirs(adata_path, exist_ok=True)
    adata.write(os.path.join(adata_path, "{sample}_RCTD_umap.h5ad"))

    sc.pl.umap(adata, color='RCTD_first_type_rat', palette=label_to_hex, show=True)
    plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
    plt.savefig(os.path.join(output_dir, f"celltype_subclass_{sample}_UMAP.png"),
                dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    