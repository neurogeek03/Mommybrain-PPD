import os
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


# Paths
project_folder = "/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists"
h5ad_dir = "/project/rrg-shreejoy/MommyBrain/Slide_tags/Pipeline_data/objects/post_bender"
metadata_path = "/scratch/mfafouti/ABC_atlas_celltyping/scratch/DW_baseline_csv_mapping_output copy.csv"
bc3_path = os.path.join(h5ad_dir, "converted_ann_data_BC3.h5ad")

print("Loading metadata...")
# Load the cell type mapping results
metadata_df = pd.read_csv(metadata_path)
print(f"Metadata shape: {metadata_df.shape}")
print(f"Metadata columns: {metadata_df.columns.tolist()}")

print("Loading BC3 data...")
# Load the BC3 h5ad file
bc3_adata = ad.read_h5ad(bc3_path)
print(f"BC3 data shape: {bc3_adata.shape}")
print(f"BC3 obs columns: {bc3_adata.obs.columns.tolist()}")

# Check if we have cell IDs that match
print(f"Metadata cell IDs: {metadata_df.iloc[:5, 0].tolist()}")
print(f"BC3 cell IDs: {bc3_adata.obs.index[:5].tolist()}")

# Merge metadata with BC3 obs
print("Merging metadata...")
# Assuming the first column is cell_id, adjust if needed
cell_id_col = metadata_df.columns[0]  # First column should be cell_id
metadata_df = metadata_df.set_index(cell_id_col)

# Merge the metadata
bc3_adata.obs = bc3_adata.obs.merge(metadata_df, left_index=True, right_index=True, how='left')

print(f"After merge - BC3 obs columns: {bc3_adata.obs.columns.tolist()}")
print(f"Number of cells with metadata: {bc3_adata.obs['subclass_name'].notna().sum()}")

# Check if subclass_name column exists
if 'subclass_name' in bc3_adata.obs.columns:
    print("Subclass names found:")
    print(bc3_adata.obs['subclass_name'].value_counts())
else:
    print("Available columns:")
    print(bc3_adata.obs.columns.tolist())
    # Try to find similar column names
    subclass_cols = [col for col in bc3_adata.obs.columns if 'subclass' in col.lower()]
    print(f"Subclass-related columns: {subclass_cols}")

print('starting computations')
sc.pp.highly_variable_genes(bc3_adata, flavor='seurat_v3', n_top_genes=2000)
sc.pp.normalize_total(bc3_adata, target_sum=1e4)
sc.pp.log1p(bc3_adata)

print('done with log normalization')
sc.pp.scale(bc3_adata, max_value=10)
sc.tl.pca(bc3_adata, svd_solver='arpack')
sc.pp.neighbors(bc3_adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(bc3_adata)
print('done with umap!')

# Create the plot
print("Creating UMAP plot...")
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# Plot UMAP colored by subclass_name
if 'subclass_name' in bc3_adata.obs.columns:
    sc.pl.umap(bc3_adata, color='subclass_name', ax=ax, show=False, 
               legend_loc='right margin', legend_fontsize=8)
    ax.set_title('BC3 Cells Colored by Subclass Name')
else:
    # If subclass_name not found, plot with a different color
    sc.pl.umap(bc3_adata, color='total_counts', ax=ax, show=False)
    ax.set_title('BC3 Cells Colored by Total Counts (subclass_name not found)')

plt.tight_layout()

# Save the plot
output_path = os.path.join(project_folder, "BC3_umap_subclass.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {output_path}")

# Save the updated h5ad file with metadata
output_h5ad = os.path.join(project_folder, "BC3_with_celltype_metadata.h5ad")
bc3_adata.write(output_h5ad)
print(f"Updated h5ad file saved to: {output_h5ad}")

plt.show()

print("Done!")
