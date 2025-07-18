import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata

# ========== PATHS ==========
in_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Spatial/anndata_objects'
out_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Integration'

metadata = os.path.join(out_dir, 'slide_seq_metadata.csv')

# ========== READING IN METADATA ==========
metadata_df = pd.read_csv(metadata)
print(metadata_df.head())

# ========== READING IN ANNDATA OBJECTS ==========
adata_list = []
for filename in os.listdir(in_dir):
    filepath = os.path.join(in_dir, filename)
    if os.path.isfile(filepath): 
        print(f"Found file: {filepath}")

    sample = filename.split('_')[0]

    # Reading in file
    ad = sc.read_h5ad(filepath)
    ad.obs['sample'] = sample
    
    # Add metadata from metadata_df
    ad.obs = ad.obs.merge(metadata_df, on='sample', how='left')

    # Filtering singlets only - OPTIONAL 
    ad = ad[ad.obs["RCTD_spot_class"] == "singlet"].copy()

    adata_list.append(ad)

    print(f'Added {filename} to the list!')

# ========== CONCATENATING ANNDATA OBJECTS ==========
adata_all = anndata.concat(
    adata_list, 
    join='inner',  # Use only shared genes
    merge='same',  # Keep .obs/.var columns only if all match
    index_unique=None
)

print(adata_all.obs.head())

adata_all.write('all_slide_seq_singlets_15.h5ad')

