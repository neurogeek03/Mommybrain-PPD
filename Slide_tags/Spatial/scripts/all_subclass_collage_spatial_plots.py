
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import pandas as pd
from pathlib import Path

# ========== PARAMETERS ==========
sample_list = ["BC13", "BC14", "BC28", "BC15" ,"BC3", "BC9"]

celltype_col = 'class_name'

# ========== PATHS ==========
mommybrain_folder = Path('/project/rrg-shreejoy/MommyBrain/Slide_tags/Pipeline_data/')
project_folder = Path('/scratch/mfafouti/Mommybrain/Slide_tags')
coords_dir = mommybrain_folder / 'spatial_coordinates'
out_dir = project_folder / 'Spatial'/'figures'
adata_dir = project_folder / 'Filtering' / 'NEW_list_merged_filtered' 
adata_path = adata_dir / 'umap_filtered_150725_slide_tags.h5ad'

fig_path = out_dir / f"{celltype_col}_NEW_combined_spatial_all_samples.png"

# ========== READING INTEGRATED ADATA FILE ==========
adata_all = sc.read_h5ad(adata_path)
print(adata_all.obs.head())
print(adata_all)

# ============ GET COLORS ============
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])
# color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)

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

# ================== LOOPING THROUGH SAMPLE LIST ==================
# Create 2x3 subplot grid
fig, axes = plt.subplots(2, 3, figsize=(20,10))

for i, sample in enumerate(sample_list): 
    # Filtering adata for the UMIs of the sample of interest
    sample_filtered_ad = adata_all[adata_all.obs['sample'] == sample].copy()
    print(sample_filtered_ad)

    group_col = celltype_col
    groups = sample_filtered_ad.obs[group_col].unique()

    df = sample_filtered_ad.obs[[celltype_col]].reset_index(names='barcode').copy()
    # df = sample_filtered_ad.obs[['barcode', celltype_col]].copy()
    df['barcode'] = df['barcode'].str.replace('-1', '', regex=False)
    print(df.head())
    print(df.shape)

    coords_path = coords_dir / f'coords_{sample}.csv'
    obs_df = pd.read_csv(coords_path)

    # Left join of spatial info 
    merged_df = df.merge(obs_df, left_on='barcode', right_on='cell_bc', how='left')
    merged_df.set_index('barcode', inplace=True)

    print(f'Shape after filtering{merged_df.shape}')
    print(merged_df.head())

    valid_obs = merged_df[(merged_df['x_um'] != 0) | (merged_df['y_um'] != 0)]
    print(valid_obs.shape)
    valid_coords = valid_obs[['x_um', 'y_um']]
    print(valid_coords.shape)
    valid_cell_subclass = valid_obs[celltype_col]

    print(valid_coords.head())
    print(valid_cell_subclass.head())

    # Plot into subplot
    row = i // 3
    col = i % 3
    ax = axes[row, col]

    sns.scatterplot(
        ax=ax,
        x=valid_coords['x_um'],
        y=valid_coords['y_um'],
        hue=valid_cell_subclass,
        palette=label_to_hex,
        s=5,
        alpha=0.7,
        legend=False  # Suppress legends in subplots
    )

    ax.set_title(f'{sample} ({len(valid_coords)} cells)', fontsize=14)
    ax.set_xlabel('X Coordinate (µm)')
    ax.set_ylabel('Y Coordinate (µm)')

# ================== PL0TTING ALL TOGETHER ==================

# Final layout and save
plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave room for the legend
plt.savefig(fig_path, dpi=300)


