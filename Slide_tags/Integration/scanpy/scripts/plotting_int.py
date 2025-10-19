
# ========== IMPORTS ==========
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import pandas as pd
import numpy as np

# ========== PATHS ==========
ad_folder = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender'
out_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Integration/scanpy'
csv_path = os.path.join(ad_folder, "hex_codes_subclass.csv")

# ========== COLOR DICTIONARY ==========
df = pd.read_csv(csv_path, header=None, names=['cell_type', 'color'])
color_dict = dict(zip(df['cell_type'], df['color']))

# ========== LOAD BBKNN ==========
h5ad_path = os.path.join(out_dir, "integrated_bbknn.h5ad")
adata_bbknn = sc.read_h5ad(h5ad_path)

# Add treatment column
adata_bbknn.obs['treatment'] = adata_bbknn.obs['sample'].apply(
    lambda x: 'OIL' if x in ['BC13', 'BC14', 'BC28'] else 'CORT'
)

# Save updated BBKNN AnnData
bbknn_out = os.path.join(out_dir, "integrated_BBKNN_t.h5ad")
adata_bbknn.write(bbknn_out)

# ========== LOAD HARMONY ==========
h5ad_path = os.path.join(out_dir, "integrated_harmony.h5ad")
adata_harmony = sc.read_h5ad(h5ad_path)

# Add treatment column
adata_harmony.obs['treatment'] = adata_harmony.obs['sample'].apply(
    lambda x: 'OIL' if x in ['BC13', 'BC14', 'BC28'] else 'CORT'
)

# Save updated Harmony AnnData
harmony_out = os.path.join(out_dir, "integrated_Harmony_t.h5ad")
adata_harmony.write(harmony_out)

# ========== PLOTTING: 2x3 UMAP GRID ==========
fig, axs = plt.subplots(2, 3, figsize=(35, 10))

# Row 0: BBKNN
sc.pl.umap(adata_bbknn, color='MapMyCells_cell_type', palette=color_dict, ax=axs[0, 0], show=False)
axs[0, 0].set_title("BBKNN - MapMyCells", fontsize=14)

sc.pl.umap(adata_bbknn, color='sample', ax=axs[0, 1], show=False)
axs[0, 1].set_title("BBKNN - Sample ID", fontsize=14)

sc.pl.umap(adata_bbknn, color='treatment', ax=axs[0, 2], show=False)
axs[0, 2].set_title("BBKNN - Treatment", fontsize=14)

# Row 1: Harmony
sc.pl.umap(adata_harmony, color='MapMyCells_cell_type', palette=color_dict, ax=axs[1, 0], show=False)
axs[1, 0].set_title("Harmony - MapMyCells", fontsize=14)

sc.pl.umap(adata_harmony, color='sample', ax=axs[1, 1], show=False)
axs[1, 1].set_title("Harmony - Sample ID", fontsize=14)

sc.pl.umap(adata_harmony, color='treatment', ax=axs[1, 2], show=False)
axs[1, 2].set_title("Harmony - Treatment", fontsize=14)

# Layout and save
plt.tight_layout()
matrix_fig_path = os.path.join(out_dir, "UMAP_summary_matrix.png")
plt.savefig(matrix_fig_path, dpi=300)