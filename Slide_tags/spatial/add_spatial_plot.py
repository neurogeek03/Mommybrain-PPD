"""
Title:        Spatial information & plot  
Description:  Adding spatial coordinates from the CurioTrekker pipileine to the .h5ad file for each sample 
Author:       Maria Eleni Fafouti 
Date:         06-05-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os
import seaborn as sns

# ========== DEFINING ESSENTIAL PATHS ==========
coords_dir = '/project/s/shreejoy/hudsonhu/RatAllCTOuts/spatial_coords_bender'
ad_folder = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender'
out_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Spatial'
csv_path = os.path.join(ad_folder, "hex_codes_subclass.csv")

df = pd.read_csv(csv_path, header=None, names=['cell_type', 'color'])
color_dict = dict(zip(df['cell_type'], df['color']))

sample_list = ['BC28', 'BC3', 'BC9', 'BC15', 'BC14', 'BC13']
# sample_list = ['BC28']

# ========== READING COORDINATE DATA AND ADDING IT AS CELL METADATA ==========
for sample in sample_list: 
    print(f"Processing {sample}...")

    # Reading coordinates file 
    coords_file = os.path.join(coords_dir, f'coords_{sample}.txt')
    coords = pd.read_csv(coords_file, header = 0, delim_whitespace=True, quotechar='"')
    print(coords.head())
    print(coords.shape)
    output_csv_file = os.path.join(out_dir, f'coords_{sample}.csv')  # Define output file path
    coords.to_csv(output_csv_file, index=False)
    # Subsetting coordinates & Barcodes only 
    coords_subset = coords.iloc[:, :3] 
    print(coords_subset.head()) 
    # Setting 1st column as index 
    coords_subset = coords_subset.set_index(coords_subset.columns[0])  # Set the first column as the index
    print(coords_subset.head()) 

    # Reading anndata file 
    sample_folder = os.path.join(ad_folder, sample)
    h5ad_path = os.path.join(sample_folder, f'_{sample}_raw_and_annotated.h5ad')
    ad = sc.read_h5ad(h5ad_path)
    print(ad)
    ad.obs.head()
    print(ad.raw.X.shape)  
    print(ad.raw.var.head()) 

    # Removing the "-1" from obs_names in the anndata obj
    ad.obs_names = ad.obs_names.str.replace('-1$', '', regex=True)

    print("Coords barcodes:", list(coords_subset.index)[:10])
    print("AnnData barcodes:", list(ad.obs_names)[:10])

    # Checking if the UMI rows match - so we get a perfect join 
    print(set(coords_subset.index).issubset(set(ad.obs_names)))

    missing = set(ad.obs_names) - set(coords_subset.index)
    print(f"Missing {len(missing)} cell IDs:", list(missing)[:10])

    # Ignoring the missing ones by adding 0, 0 to the spatial coords 
    ad.obs['x_um'] = coords_subset['x_um'].reindex(ad.obs_names, fill_value=0)
    ad.obs['y_um'] = coords_subset['y_um'].reindex(ad.obs_names, fill_value=0)

    # Preview the entire obs DataFrame, including the new columns
    pd.set_option('display.max_columns', None)
    print(ad.obs.head())
    obs_df = ad.obs.copy()
    obs_df['cell_bc'] = ad.obs_names

    # Saving anndata obj
    new_file_dir = os.path.join(ad_folder, f"{sample}", f"{sample}_spatial.h5ad")
    ad.write(new_file_dir)

     # ========== EXTRACTING SPATIAL DATA ==========
    # Fitering out barcodes with no spatial info
    valid_obs = obs_df[(obs_df['x_um'] != 0) | (obs_df['y_um'] != 0)]

    valid_coords = valid_obs[['x_um', 'y_um']]  # Extract the spatial coordinates
    valid_cell_class = valid_obs['MapMyCells_cell_type']  # Extract the cell class
    valid_cell_bcs = valid_obs['cell_bc']  # Extract the cell barcodes

    # Ensuring the dfs are of the same length
    assert len(valid_coords) == len(valid_cell_class)
    assert len(valid_coords) == len(valid_cell_bcs) 

    print(f'We will plot: {len(valid_coords)} coordinate pairs!')

    # ========== CREATING SPATIAL PLOT ==========
    figure_path = os.path.join(out_dir, f'{sample}_spatial_cell_class.png')
    plt.figure(figsize=(8, 8))

    # Plot the spatial locations, coloring by cell class
    sns.scatterplot(
        x=valid_coords['x_um'], 
        y=valid_coords['y_um'], 
        hue=valid_cell_class, 
        palette=color_dict, 
        s=5, 
        alpha=0.7)

    # Style
    plt.xlabel('X Coordinate (um)')
    plt.ylabel('Y Coordinate (um)')
    plt.title(f'{sample} - Spatial Distribution of {len(valid_coords)} coordinate pairs')
    plt.legend(title='Cell Class', bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)
    plt.savefig(figure_path, dpi=300, bbox_inches='tight', pad_inches=0.1)



    



