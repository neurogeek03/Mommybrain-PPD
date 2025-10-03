import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata

# ========== PARAMS ==========
singlet_score_trheshold = 0

# ========== PATHS ==========
in_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/OUT/FINAL_RCTD_newgenelist/anndata_objects' #CHANGE
out_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist' #CHANGE
metadata = '/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR/slide_seq_metadata.csv'
os.makedirs(out_dir, exist_ok=True)

# ========== READING IN METADATA ==========
metadata_df = pd.read_csv(metadata)
print(metadata_df.head())

# ========== READING IN ANNDATA OBJECTS ==========
adata_list = []
for filename in os.listdir(in_dir):
    filepath = os.path.join(in_dir, filename)
    if os.path.isfile(filepath): 
        print(f"Found file: {filepath}")

    sample = filename.split("_")[0]
    print(f'Processing sample {sample}..')

    # Reading in file
    ad = sc.read_h5ad(filepath)
    ad.obs['sample'] = sample
    
    # Add metadata from metadata_df
    ad.obs = ad.obs.merge(metadata_df, on='sample', how='left')

    # Filtering singlets only - OPTIONAL 
    ad = ad[ad.obs["RCTD_spot_class_rat"] == "singlet"].copy()

    # Filtering singlets score = 330
    ad = ad[ad.obs["RCTD_singlet_score_rat"] > singlet_score_trheshold].copy()
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

adata_all.obs_names_make_unique() # some barcodes are shared across samples 

print(adata_all.obs.dtypes)
print(adata_all.obs["pregnancy"].unique())
print(adata_all.obs["pregnancy"].map(type).value_counts())


for col in ["pregnancy", "day", "treatment", "sample"]:
    adata_all.obs[col] = adata_all.obs[col].astype(str)


adata_all.write(os.path.join(out_dir,f'NEW_genelist_singlet_score_{singlet_score_trheshold}_slide_seq_15.h5ad'))

