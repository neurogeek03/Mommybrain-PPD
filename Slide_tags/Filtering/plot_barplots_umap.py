"""
Title: Plotting results of filtered data
Description:  Umap & barplots (at classs & subclass level)
Author:   Maria Eleni Fafouti 
Date: 12-05-2025
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
out_dir = os.path.join(project_path,'Filtering')
adata_file = os.path.join(out_dir, "FILTERED_persample_subclass_doublets_harmony.h5ad")

n_cells = 50

# ========== READING IN DATA  ==========
adata = sc.read_h5ad(adata_file)

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

# ========== PLOTTING ==========

adata.obs['MapMyCells_subclass_name'] = pd.Categorical(
    adata.obs['MapMyCells_subclass_name'],
    categories=my_order,
    ordered=True
)

# Update color info in AnnData
adata.uns['MapMyCells_subclass_name_colors'] = [
    color_dict.get(subcls, '#cccccc') for subcls in my_order
]

subclass_counts = adata.obs['MapMyCells_subclass_name'].value_counts().reindex(my_order)

# -------------------------
# Assign colors for subclasses in AnnData
# -------------------------
# subclass_categories = adata.obs['MapMyCells_subclass_name'].cat.categories
# subclass_colors = [color_dict.get(subcls, '#cccccc') for subcls in subclass_categories]
# adata.uns['MapMyCells_subclass_name_colors'] = subclass_colors

# -------------------------
# Plot: UMAP colored by subclass
# -------------------------
sc.set_figure_params(figsize=(8, 6))

# Create a figure with 1 row and 2 columns
fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Adjust width as needed

custom_palette = ['#1f77b4', '#ff7f0e']

# Second UMAP: colored by treatment
sc.pl.umap(
    adata,
    color='treatment',
    ax=axes[0],
    show=False,
    legend_loc=None,
    palette=custom_palette,
    title="Treatment"
)

# First UMAP: colored by subclass name
sc.pl.umap(
    adata,
    color='MapMyCells_subclass_name',
    ax=axes[1],
    show=False,
    legend_loc=None,
    title=f"Subclass name (n = {adata.n_obs})"
)

# Adjust layout and save
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'side_by_side_umaps.png'), dpi=300)

# -------------------------
# Plot: number of cells per subclass
# -------------------------
subclass_counts_sorted = subclass_counts.sort_values(ascending=True)

# categories = adata.obs['MapMyCells_subclass_name'].cat.categories
# palette = adata.uns['MapMyCells_subclass_name_colors']

# # Map categories to colors
# color_map = dict(zip(categories, palette))

# # Sort subclass counts however you want, then get colors for those subclasses
# subclass_counts_sorted = subclass_counts.sort_values(ascending=True)
# bar_colors = [color_map.get(subcls, '#cccccc') for subcls in subclass_counts_sorted.index]

# # Plot
# plt.figure(figsize=(10, 12))
# plt.barh(subclass_counts_sorted.index, subclass_counts_sorted.values, color=bar_colors)
# plt.xlabel('Number of Cells')
# plt.title(f'Cell Number per Subclass (n ≥ {n_cells})')
# plt.grid(False)  # turns off the grid
# plt.gca().invert_yaxis()
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, 'barplot_cells_per_subclass.png'), dpi=300)

plt.figure(figsize=(10, 12))
plt.barh(
    subclass_counts_sorted.index, 
    subclass_counts_sorted.values,
    color=[color_dict.get(subcls, '#cccccc') for subcls in subclass_counts_sorted.index]
)
plt.xlabel('Number of Cells')
plt.title(f'Cell Number per Subclass (n ≥ {n_cells})')
plt.grid(False)  # turns off the grid

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'new_barplot_cells_per_subclass.png'), dpi=300)
plt.close()


