"""
Title: Doublet Visualization
Description:  Visualizing doublets in the UMAP space and after running scDblFinder 
Author:   Maria Eleni Fafouti 
Date: 09-05-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import anndata
import scanpy as sc
import pandas as pd
import numpy as np

# ========== DEFINING ESSENTIAL PATHS ==========
project_path = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags'
ad_file = os.path.join(project_path,'Integration/scanpy/integrated_Harmony_t.h5ad')
out_dir = os.path.join(project_path,'Doublet_detection/scanpy')
csv_folder = os.path.join(project_path, 'Doublet_detection/scanpy/tmp_dir')

adata = sc.read_h5ad(ad_file)

# Make a copy of obs with a composite index
adata.obs["merge_key"] = adata.obs["sample"].astype(str) + "_" + adata.obs_names.astype(str)

print(adata.obs_names[:5]) #'TGTGCCTCAACGTTCC-1' # FIXME

samples = adata.obs["sample"].unique()

print(adata)

all_annotations = []

for sample in samples:
    csv_path = os.path.join(csv_folder, f'{sample}_scdblfinder.csv')
    dbl = pd.read_csv(csv_path)
    
    # Ensure the barcode column is named "barcode"
    dbl = dbl.rename(columns={dbl.columns[0]: "barcode"})
    
    # Add a 'sample' column to match adata.obs
    dbl["sample"] = sample

    dbl["merge_key"] = dbl["sample"].astype(str) + "_" + dbl["barcode"].astype(str)

    all_annotations.append(dbl)

# Combine all sample annotations into one DataFrame
annotations_df = pd.concat(all_annotations)

print(annotations_df.head(20))

# Merge with adata.obs
adata.obs = adata.obs.merge(annotations_df, on="merge_key", how="left")

# (Optional) remove the merge key if no longer needed
adata.obs.drop(columns=["merge_key"], inplace=True)

# ========== SAVING ANNDATA AS H5AD ==========
new_file_dir = os.path.join(out_dir, "doublets_harmony.h5ad")
adata.write(new_file_dir)

# ========== VISUALZING DOUBLETS IN UMAP SPACE ==========
sc.pl.umap(adata, color='scDblFinder.class', palette=['#1f77b4', '#d62728'], title='Doublet Detection by scDblFinder')
fig_path = os.path.join(out_dir, "UMAP_doublets2.png")
plt.tight_layout()
plt.savefig(fig_path, dpi=300)

# ========== BARPLOT OF DOUBLETS/SINGLETS PER ALLEN INSTITUTE CELL CLASS ==========
df_all = adata.obs[[
    'scDblFinder.class',
    'MapMyCells_cell_type'
]].copy()
df_all['obs_names'] = adata.obs_names

# Count total cells per type
total_counts = df_all['MapMyCells_cell_type'].value_counts().rename('total_cells')

# Count doublets per type
doublet_counts = df_all[df_all['scDblFinder.class'] == 'doublet'] \
                    .groupby('MapMyCells_cell_type') \
                    .size().rename('doublet_cells')

# Combine counts into one DataFrame
counts_df = pd.concat([total_counts, doublet_counts], axis=1).fillna(0)
counts_df['singlet_cells'] = counts_df['total_cells'] - counts_df['doublet_cells']

# Sorting
counts_df = counts_df.sort_values('total_cells', ascending=False)

# Optional: log-transform
counts_df['log_total'] = np.log1p(counts_df['total_cells'])
counts_df['log_doublet'] = np.log1p(counts_df['doublet_cells'])
counts_df['log_singlet'] = np.log1p(counts_df['singlet_cells'])

# Plot stacked bar
plt.figure(figsize=(10, 8))

cell_types = counts_df.index
singlet_vals = counts_df['singlet_cells']
doublet_vals = counts_df['doublet_cells']

plt.barh(cell_types, singlet_vals, color='gray', label='Singlets')
plt.barh(cell_types, doublet_vals, left=singlet_vals, color='red', label='Doublets')

plt.xlabel('Number of Cells')
plt.title('Doublet Proportion per Cell Type')
plt.legend()
plt.tight_layout()

fig_path = os.path.join(out_dir, "stacked_barplot_doublet_proportions.png")
plt.savefig(fig_path, dpi=300)
plt.close()
