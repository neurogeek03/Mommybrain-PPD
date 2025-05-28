"""
Title: Cell type filtering & plotting
Description:  Filetring cells based on a threshold for subclass and plotting in the UMAP space
Author:   Maria Eleni Fafouti 
Date: 09-05-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm  # for colormaps
from matplotlib.colors import Normalize
import os
import anndata as ad
import pandas as pd
import numpy as np

# ========== DEFINING ESSENTIAL PATHS ==========
project_path = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags'
adata_dir = os.path.join(project_path, "Post_bender") # path with adata obj for individual samples
ad_file = os.path.join(project_path,'Doublet_detection/scanpy/doublets_harmony.h5ad')
out_dir = os.path.join(project_path,'Filtering')

# ========== READING IN DATA  ==========
adata = sc.read_h5ad(ad_file)
print(adata.obs.head())
print(adata)

# Refining obs columns 
del adata.obs['MapMyCells_broad_type']

# Renaming
adata.obs.rename(columns={'MapMyCells_cell_type':'MapMyCells_class_name' }, inplace=True)
adata.obs.rename(columns={'barcode_x':'barcode' }, inplace=True)
adata.obs.rename(columns={'sample_x':'sample' }, inplace=True)

print(adata.obs[['sample', 'barcode']].head())

# ========== ADDING VAR INFO ==========
var_adata_path = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender/BC9/converted_ann_data_BC9.h5ad"

adata_older = sc.read_h5ad(var_adata_path)
# print(adata_older.var.head())

adata.var = adata_older.var.copy()
# print(adata.var.head())

#TODO join subclass based on sample and barcode combination. Each sample has a separate csv showing which subclass each cluster belongs to

# ========== ADDING SUBCLASS INFO TO OBS ==========

sample_list = [ "BC14", "BC28", "BC13", "BC15", "BC3", "BC9"]

# Create an empty list to store dataframes
annot_dfs = []

for sample in sample_list:
    print(f'Sample {sample} is being processed')

    # Load the MapMyCells annotations
    csv_path = os.path.join(adata_dir, f"{sample}", f'{sample}_mapmycells.csv')
    annotations = pd.read_csv(csv_path, skiprows=0, index_col=0)  # adjust index_col if needed

    # Check the first few rows
    print(annotations.head())
    # print(annotations.index[:5])

    annotations['sample'] = sample
    annotations['barcode'] = annotations.index

    # Select relevant columns
    annot_dfs.append(annotations[['sample', 'barcode', 'subclass_name', 'subclass_bootstrapping_probability']])

all_annotations = pd.concat(annot_dfs)

# Merge with adata.obs using sample and barcode
adata.obs = adata.obs.merge(all_annotations,
                            on=['sample', 'barcode'],
                            how='left')

print(adata)

print(adata.obs.head())

# Optional: rename the column for clarity
adata.obs.rename(columns={'subclass_name': 'MapMyCells_subclass_name'}, inplace=True)
adata.obs.rename(columns={'subclass_bootstrapping_probability': 'MapMyCells_subclass_bootstrapping_probability'}, inplace=True)

print(adata.obs.head())

# # Saving object
# out_path = os.path.join(out_dir, "non_filtered_subclass_doublets_harmony.h5ad")
# adata.write(out_path)

# ========== FILTERING OUT DOUBLETS ==========
adata_filtered = adata[adata.obs['scDblFinder.class'] == "singlet"].copy()


# # ========== FILTERING BASED ON TOTAL CELLS ==========
# # Filtering number of cells per subclass 
# subclass_counts_all = adata_filtered.obs['MapMyCells_subclass_name'].value_counts()

# # Keep only subclasses with n cells
# n_cells = 300
# valid_subclasses = subclass_counts_all[subclass_counts_all >= n_cells].index

# # Filter adata accordingly
# adata_filtered = adata_filtered[adata_filtered.obs['MapMyCells_subclass_name'].isin(valid_subclasses)].copy()

# ========== FILTERING PER SAMPLE ==========
# Compute number of cells per (sample, subclass) group
group_counts = (
    adata_filtered.obs
    .groupby(['MapMyCells_subclass_name', 'sample'])
    .size()
    .unstack(fill_value=0)  # rows: subclass, columns: sample
)

# Keep only subclasses with ≥50 cells in *all* samples
n_cells = 50
valid_subclasses = group_counts[(group_counts >= n_cells).all(axis=1)].index

adata_filtered = adata_filtered[
    adata_filtered.obs['MapMyCells_subclass_name'].isin(valid_subclasses)
    ].copy()

out_path = os.path.join(out_dir, "FILTERED_persample_subclass_doublets_harmony.h5ad")
adata_filtered.write(out_path)

# ========== PLOTTING ==========
# Allen hex codes for classes
csv_path = os.path.join(project_path, "Post_bender", "hex_codes_subclass.csv")
df = pd.read_csv(csv_path, header=None, names=['cell_type', 'color'])
color_dict = dict(zip(df['cell_type'], df['color']))

# cells per class/subclass 
# Create a copy of relevant columns from adata.obs
df_all = adata_filtered.obs[['MapMyCells_class_name', 'MapMyCells_subclass_name']].copy()
df_all['obs_names'] = adata_filtered.obs_names

# Count cells per subclass
subclass_counts_all = df_all['MapMyCells_subclass_name'].value_counts()

# Recalculate counts
subclass_counts = df_all['MapMyCells_subclass_name'].value_counts().sort_values(ascending=True)
class_counts = df_all['MapMyCells_class_name'].value_counts().sort_values(ascending=True)

class_colors = [color_dict.get(cls, '#cccccc') for cls in class_counts.index]
# -------------------------
# Plot: number of cells per class
# -------------------------
plt.figure(figsize=(10, 6))
plt.barh(class_counts.index, class_counts.values, color=class_colors)
plt.xlabel('Number of Cells')
plt.title(f'Cell Number per Class (≥{n_cells} cells per subclass)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'barplot_cells_per_class.png'), dpi=300)
plt.close()

subclass_to_class = df_all.set_index('MapMyCells_subclass_name')['MapMyCells_class_name'].to_dict()
subclass_colors = [
    color_dict.get(subclass_to_class.get(subcls, ''), '#cccccc')
    for subcls in subclass_counts.index
]

# -------------------------
# Plot: number of cells per subclass
# -------------------------
plt.figure(figsize=(10, 12))  # More height for subclasses
plt.barh(subclass_counts.index, subclass_counts.values, color=subclass_colors)
plt.xlabel('Number of Cells')
plt.title(f'Cell Number per Subclass (n ≥ {n_cells})')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'barplot_cells_per_subclass.png'), dpi=300)
plt.close()

sc.set_figure_params(figsize=(12,6))
# adata_filtered.uns['MapMyCells_subclass_name_colors'] = sns.color_palette("Set3", as_cmap=False).as_hex()
sc.pl.umap(
    adata_filtered,
    color='MapMyCells_subclass_name',
    show=False
)

# Add title and save
plt.gca().set_title(f"UMAP - MapMyCells class name (n = {adata_filtered.n_obs})", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'filtered_umap_MapMyCells_Harmony.png'), dpi=300)
plt.close()
