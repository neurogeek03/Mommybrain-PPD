#!/usr/bin/env python

import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse
import os
# import imageio
from matplotlib.patches import Patch

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

# ========== ARGS ==========
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help="Path to input annotated anndata objects w/ RCTD output", default ='/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/OUT/B03_default_full_genes/anndata_objects')
parser.add_argument("-o", "--output_dir", help="Path to output folder for the figures", default='/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/Filtering/OUT/gif')
args = parser.parse_args()

# ========== PATHS ==========
adata_dir = args.input_dir
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

# ========== GET FILES ==========
h5ad_files = sorted([f for f in os.listdir(adata_dir) if f.endswith(".h5ad")])
spot_class_column = "RCTD_spot_class_rat"
type_column = "RCTD_first_type_rat"
filter_column = "RCTD_singlet_score_rat"
filter_values = list[0, 100]

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

images = []

# ========== LOOP OVER FILES ==========
for h5ad_file in h5ad_files:    
    for filter_value in filter_values:
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

        adata = adata[adata.obs[filter_column]>filter_value].copy()
        num_cells = adata.n_obs

        print("➡️ Starting computations...")
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, zero_center=False)
        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        sc.tl.umap(adata, random_state=42)
        print(f"✅ Done with UMAP for {sample}")

        # Count cell types
        type_counts = adata.obs[type_column].value_counts()
        top_types = type_counts.head(30).index.tolist()

        # Build legend elements using your label_to_hex
        legend_elements = [
            Patch(facecolor=label_to_hex.get(t, "#808080"), label=t) for t in top_types
        ]

        sc.pl.umap(adata, color=type_column, palette=label_to_hex, show=True)
        plt.title(f"UMAP of singlets (n = {num_cells}) after RCTD mouse", fontsize=14)
        # Build custom legend for top 30 only
        legend_elements = [
            Patch(facecolor=label_to_hex.get(t, "#808080"), label=t) for t in top_types
        ]
        plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5), title="Top 30 subclasses", fontsize=7)
        img_path = os.path.join(output_dir, f"min_{filter_value}_{sample}_UMAP.png")
        plt.savefig(img_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()
        images.append(img_path)
        

# gif_path = os.path.join(output_dir, "umap_animation.gif")
# with imageio.get_writer(gif_path, mode='I', duration=3) as writer:
#     for img in images:
#         frame = imageio.imread(img)
#         writer.append_data(frame)

# print(f"✅ GIF saved at {gif_path}")
