"""
Title:        Spatial Plots colored by Allen Brain Institute Subclass
Description:  In the integrated data file: search for each sample, add spatial coords and plot. 
Author:       Maria Eleni Fafouti 
Date:         28-05-2025
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


# ========== PATHS ==========
coords_dir = '/project/s/shreejoy/hudsonhu/RatAllCTOuts/spatial_coords_bender'
out_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Spatial'
project_path = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags"
adata_dir = os.path.join(project_path, "Post_bender")
filtered_ad = os.path.join(project_path, "Filtering", "non_filtered_subclass_doublets_harmony.h5ad") # can change this to filtered adata as well 

# ========== READING INTEGRATED ADATA FILE ==========
adata_all = sc.read_h5ad(filtered_ad)
print(adata_all.obs.head())
print(adata_all)

# ========== PARAMETERS ==========
sample_list = ["BC13", "BC14", "BC28", "BC15" ,"BC3", "BC9"]

color_dict = {
    "001 CLA-EPd-CTX Car3 Glut": "#ffff00",        # yellow
    "004 L6 IT CTX Glut": "#00ffff",                # cyan
    "006 L4/5 IT CTX Glut": "#ff00ff",              # magenta
    "007 L2/3 IT CTX Glut": "#ff0000",              # red
    "009 L2/3 IT PIR-ENTl Glut": "#006666",         # teal-ish dark cyan
    "012 MEA Slc17a7 Glut": "#0080ff",               # bright blue
    "014 LA-BLA-BMA-PA Glut": "#660000",             # dark maroon
    "016 CA1-ProS Glut": "#33cc33",                  # bright green
    "017 CA3 Glut": "#663300",                        # brown
    "022 L5 ET CTX Glut": "#000099",                  # navy blue
    "025 CA2-FC-IG Glut": "#00ffcc",                   # bright turquoise
    "030 L6 CT CTX Glut": "#cc9966",                  # tan/light brown
    "037 DG Glut": "#003366",                          # dark blue
    "045 OB-STR-CTX Inh IMN": "#9999ff",              # lavender/light purple
    "061 STR D1 Gaba": "#666666",                      # medium gray
    "062 STR D2 Gaba": "#330000",                      # very dark red/burgundy
    "093 RT-ZI Gnb3 Gaba": "#666666",                  # gray (looks same as 061, but different tone)
    "101 ZI Pax6 Gaba": "#330000",                      # dark maroon (same as 062)
    "151 TH Prkcd Grin2c Glut": "#003300",             # dark green
    "152 RE-Xi Nox4 Glut": "#ffd9e6",                   # light pink
    "318 Astro-NT NN": "#0000ff",                       # blue
    "319 Astro-TE NN": "#666699",                       # muted purple-gray
    "325 CHOR NN": "#ff007f",                           # pink/magenta
    "326 OPC NN": "#666666",                            # gray (same as above grays)
    "327 Oligo NN": "#990000",                          # dark red
    "330 VLMC NN": "#999900",                           # olive green
    "331 Peri NN": "#33cccc",                           # teal/turquoise
    "333 Endo NN": "#ffcc99",                           # peach/light orange
    "334 Microglia NN": "#ffb6c1",                      # light pink
    "338 Lymphoid NN": "#cc33cc",                       # purple-pink
}

my_order = [
    "007 L2/3 IT CTX Glut",
    "006 L4/5 IT CTX Glut",
    "004 L6 IT CTX Glut",
    "009 L2/3 IT PIR-ENTl Glut",
    "030 L6 CT CTX Glut",
    "022 L5 ET CTX Glut",
    "001 CLA-EPd-CTX Car3 Glut",
    "014 LA-BLA-BMA-PA Glut",
    "012 MEA Slc17a7 Glut",
    "025 CA2-FC-IG Glut",
    "017 CA3 Glut",
    "016 CA1-ProS Glut",
    "037 DG Glut",
    "152 RE-Xi Nox4 Glut",
    "151 TH Prkcd Grin2c Glut",
    "061 STR D1 Gaba",
    "062 STR D2 Gaba",
    "093 RT-ZI Gnb3 Gaba",
    "101 ZI Pax6 Gaba",
    "045 OB-STR-CTX Inh IMN",
    "318 Astro-NT NN",
    "319 Astro-TE NN",
    "326 OPC NN",
    "327 Oligo NN",
    "330 VLMC NN",
    "331 Peri NN",
    "333 Endo NN",
    "334 Microglia NN",
    "338 Lymphoid NN",
    "325 CHOR NN"
]

adata_all.obs['MapMyCells_subclass_name'] = pd.Categorical(
    adata_all.obs['MapMyCells_subclass_name'],
    categories=my_order,
    ordered=True
)


# ================== LOOPING THROUGH SAMPLE LIST ==================
# Create 2x3 subplot grid
fig, axes = plt.subplots(2, 3, figsize=(20,10))

for i, sample in enumerate(sample_list): 
    # Filtering adata for the UMIs of the sample of interest
    sample_filtered_ad = adata_all[adata_all.obs['sample'] == sample].copy()
    print(sample_filtered_ad)

    group_col = "MapMyCells_subclass_name"
    groups = sample_filtered_ad.obs[group_col].unique()

    df = sample_filtered_ad.obs[['barcode', 'MapMyCells_subclass_name']].copy()
    df['barcode'] = df['barcode'].str.replace('-1', '', regex=False)
    print(df.head())
    print(df.shape)

    # Reading individual sample w/ spatial coords
    adata_file = os.path.join(adata_dir, f"{sample}", f"{sample}_spatial.h5ad")
    ad = sc.read_h5ad(adata_file)
    ad.obs['barcode'] = ad.obs_names

    print(f"Shape before filering: {ad.obs.shape}")
    print(ad.obs.head())

    # Extract spatial data from obs
    obs_df = ad.obs.copy()

    # Left join of spatial info 
    merged_df = df.merge(obs_df, on='barcode', how='left')
    merged_df.set_index('barcode', inplace=True)

    print(f'Shape after filtering{merged_df.shape}')
    print(merged_df.head())

    valid_obs = merged_df[(merged_df['x_um'] != 0) | (merged_df['y_um'] != 0)]
    print(valid_obs.shape)
    valid_coords = valid_obs[['x_um', 'y_um']]
    print(valid_coords.shape)
    valid_cell_subclass = valid_obs['MapMyCells_subclass_name']

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
        palette=color_dict,
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
plt.savefig(os.path.join(out_dir, "unfiltered_combined_spatial_all_samples.png"), dpi=300)


