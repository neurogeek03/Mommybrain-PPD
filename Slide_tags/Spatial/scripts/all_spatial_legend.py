# ========== IMPORTS ==========
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import pandas as pd
from pathlib import Path
from matplotlib.patches import Patch  # Make sure this is imported

# ========== PARAMETERS ==========
sample_list = ["BC13", "BC14", "BC28", "BC15" ,"BC3", "BC9"]

celltype_col = 'subclass_name'

# ========== PATHS ==========
mommybrain_folder = Path('/project/rrg-shreejoy/MommyBrain/Slide_tags/Pipeline_data/')
project_folder = Path('/scratch/mfafouti/Mommybrain/Slide_tags')
coords_dir = mommybrain_folder / 'spatial_coordinates'
out_dir = project_folder / 'Spatial'/'figures'
adata_dir = project_folder / 'Filtering' / 'out'
adata_path = adata_dir / "PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"

# Make sure output directory exists
out_dir.mkdir(parents=True, exist_ok=True)
fig_path = out_dir / f"mt_filtered_{celltype_col}_NEW_combined_spatial_all_samples.png"

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

# ============ GET TOP 30 MOST ABUNDANT CELL TYPES ============
# Get the names of the top 30 most frequent cell types from the whole dataset
top_30_types = adata_all.obs[celltype_col].value_counts().head(30).index.tolist()
print(f"Top 30 types identified for legend: {top_30_types}")

# ================== LOOPING THROUGH SAMPLE LIST ==================
# Create 2x3 subplot grid
# --- 1. CHANGED: Increased figure height for the long legend ---
fig, axes = plt.subplots(2, 3, figsize=(20, 15)) 

for i, sample in enumerate(sample_list): 
    # Filtering adata for the UMIs of the sample of interest
    sample_filtered_ad = adata_all[adata_all.obs['sample'] == sample].copy()
    print(f"Processing {sample} with {sample_filtered_ad.n_obs} cells")

    group_col = celltype_col
    groups = sample_filtered_ad.obs[group_col].unique()

    df = sample_filtered_ad.obs[[celltype_col]].reset_index(names='barcode').copy()
    df['barcode'] = df['barcode'].str.replace('-1', '', regex=False)

    coords_path = coords_dir / f'coords_{sample}.csv'
    if not coords_path.exists():
        print(f"Warning: Coords file not found for {sample}. Skipping plot.")
        # Set title and turn off axis for the empty plot
        row = i // 3
        col = i % 3
        ax = axes[row, col]
        ax.set_title(f'{sample} (No coords file)', fontsize=14)
        ax.axis('off')
        continue # Skip this sample
        
    obs_df = pd.read_csv(coords_path)

    # Left join of spatial info 
    merged_df = df.merge(obs_df, left_on='barcode', right_on='cell_bc', how='left')
    merged_df.set_index('barcode', inplace=True)

    valid_obs = merged_df[(merged_df['x_um'] != 0) | (merged_df['y_um'] != 0)]
    valid_coords = valid_obs[['x_um', 'y_um']]
    valid_cell_subclass = valid_obs[celltype_col]

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

# ============ CREATE AND PLOT LEGEND ============

# 1. Create a list of legend "handles" (Patches) for the top 30 types
legend_patches = []
for cell_type in top_30_types:
    color = label_to_hex.get(cell_type, '#cccccc') # Use .get() for safety
    patch = Patch(color=color, label=cell_type)
    legend_patches.append(patch)

# 2. Add the custom legend to the figure
fig.legend(
    handles=legend_patches,
    # --- 2. CHANGED: Positioned legend anchor at the top-left ---
    loc='upper left',             
    # --- 3. CHANGED: Placed anchor near top-right of figure box ---
    bbox_to_anchor=(0.85, 0.95),  
    ncols=1,                     # --- 4. CHANGED: Set to one column ---
    title=f"Top 30 {celltype_col}",
    fontsize='small',
    title_fontsize='medium'
)

# ================== PL0TTING ALL TOGETHER ==================

# Final layout and save
plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave room for the legend (this rect is still correct)
plt.savefig(fig_path, dpi=300)
print(f"Saved figure with legend to {fig_path}")
plt.close(fig) # Close the figure to free up memory