import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata

# ========== PATHS ==========
in_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/data'
out_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference_coronal'
os.makedirs(out_dir, exist_ok=True)

metadata_csv = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/subset_cellids_celltypes_metadata_20241115.csv'
min_genes = 200
min_counts = 500
target_cells_per_region = 50000
min_cells_per_type = 25
adata_list = []

# ========== LOAD METADATA ==========
print("Loading metadata...")
meta_df = pd.read_csv(metadata_csv)
meta_df["region"] = meta_df["feature_matrix_label"].str.split("-").str[2]

print(f"Metadata loaded with {meta_df.shape[0]} rows")
print(meta_df.head())

# ========== PROCESS EACH REGION ==========
for filename in os.listdir(in_dir):
    if not filename.endswith('.h5ad'):
        continue

    filepath = os.path.join(in_dir, filename)
    brain_region = filename.split('-', 2)[-1].replace("-raw.h5ad", "")
    print(f"\n🧠 Processing region: {brain_region}")

    # Load regional AnnData files from ABC atlas
    ad = sc.read_h5ad(filepath)
    print(f"Loaded {ad.n_obs} cells")
    print(ad.obs.columns)
    print(ad.obs_names[:5])
    print(meta_df['cell_label'].head())

    # Filter low-quality cells
    sc.pp.filter_cells(ad, min_genes=min_genes)
    sc.pp.filter_cells(ad, min_counts=min_counts)
    print(f"Filtered: {ad.n_obs} cells remain")

    # Downsample
    if ad.n_obs > target_cells_per_region:
        sampled_idx = np.random.choice(ad.obs_names, size=target_cells_per_region, replace=False)
        ad = ad[sampled_idx].copy()
        print(f"Downsampled to {ad.n_obs} cells")

    # Add region label
    ad.obs['region'] = brain_region

    # ========== MERGE METADATA ==========
    # Ensure ad.obs index is cell barcodes matching meta_region cell_label
    ad.obs = ad.obs.copy()
    ad.obs.index = ad.obs_names

    #Filter metadata for current region and matching barcodes
    meta_region = meta_df[
        (meta_df['cell_label'].isin(ad.obs_names)) &
        (meta_df['region'].isin([brain_region]))
        ].set_index("cell_label")
    
    pd.set_option('display.max_columns', None)  # Show all columns
    pd.set_option('display.width', 0)  # Prevent line wrapping
    print(meta_region.head(20))

    # Drop overlapping columns in ad.obs if any
    overlapping_cols = ad.obs.columns.intersection(meta_region.columns)
    if len(overlapping_cols) > 0:
        print(f"Dropping overlapping columns from ad.obs: {list(overlapping_cols)}")
        ad.obs.drop(columns=overlapping_cols, inplace=True)

    # Join metadata into ad.obs by index
    ad.obs = ad.obs.join(meta_region, how='left')

    print(f"ad.obs shape after join: {ad.obs.shape}")
    print(ad.obs.head())
    print("\n📋 Columns in adata_all.obs:")
    for col in ad.obs.columns:
        print(f"\n🧷 {col}:\n{ad.obs[col].value_counts(dropna=False)}")

    adata_list.append(ad)

# ========== CONCATENATE ==========
print("\n🔗 Concatenating all regions...")
adata_all = anndata.concat(
    adata_list,
    join='outer',
    merge='same',
    index_unique=None
)
print(f"Final concatenated shape: {adata_all.shape}")

print("\n📋 Columns in adata_all.obs:")
for col in adata_all.obs.columns:
    print(f"\n🧷 {col}:\n{adata_all.obs[col].value_counts(dropna=False)}")

# ========== FILTER CELL TYPES BY COUNT ==========
# Filter out cell types with fewer than min_cells_per_type cells
cell_type_column = 'subclass'  # <<< IMPORTANT: Check if this is the correct column name for cell types
print(f"\n🔬 Filtering cell types in '{cell_type_column}' with less than {min_cells_per_type} cells...")

# Check if the column exists
if cell_type_column not in adata_all.obs.columns:
    print(f"⚠️  Warning: Column '{cell_type_column}' not found in adata_all.obs. Skipping cell type filtering.")
    print(f"    Available columns are: {adata_all.obs.columns.tolist()}")
else:
    # Count cells per type
    cell_counts = adata_all.obs[cell_type_column].value_counts()

    # Identify cell types to keep
    valid_cell_types = cell_counts[cell_counts >= min_cells_per_type].index

    print(f"Found {len(valid_cell_types)} out of {len(cell_counts)} cell types with >= {min_cells_per_type} cells.")

    # Filter the AnnData object
    n_obs_before = adata_all.n_obs
    adata_all = adata_all[adata_all.obs[cell_type_column].isin(valid_cell_types)].copy()
    print(f"Filtered from {n_obs_before} to {adata_all.n_obs} cells.")

# ========== CREATE AND SAVE SUMMARY CSV ==========
print("\n📄 Creating and saving summary CSV...")
if cell_type_column in adata_all.obs.columns and 'region' in adata_all.obs.columns:
    # Create a dataframe with the required columns from the filtered AnnData
    summary_df = adata_all.obs[['region', cell_type_column]].copy()
    summary_df.index.name = 'cell_label'
    summary_df.reset_index(inplace=True)

    # Calculate total number of cells per subclass
    subclass_counts = summary_df[cell_type_column].value_counts()

    # Add the counts to the summary dataframe
    summary_df['total_cells_per_subclass'] = summary_df[cell_type_column].map(subclass_counts)

    # Save to CSV
    csv_output_path = os.path.join(out_dir, "cell_subclass_summary.csv")
    summary_df.to_csv(csv_output_path, index=False)
    print(f"✅ Saved summary with {len(summary_df)} cells to: {csv_output_path}")
    print("CSV preview:")
    print(summary_df.head())
else:
    print(f"⚠️  Warning: Could not create summary CSV because required columns ('{cell_type_column}', 'region') were not found.")

# ========== SAVE ==========
output_filepath = os.path.join(out_dir, "celltypes_WMB_high_qual_cells_all_regions_max50k.h5ad")
adata_all.write_h5ad(output_filepath)
print(f"✅ Saved to: {output_filepath}")
