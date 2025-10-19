import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata 
from pathlib import Path

# ========== PATHS ==========
project_folder = Path ('/scratch/mfafouti/ABC_atlas_celltyping')
in_dir = project_folder / 'celltyped_doublets'
out_dir = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Filtering/NEW_list_merged_filtered')
os.makedirs(out_dir, exist_ok=True)

# ========== READING IN ANNDATA OBJECTS ==========
adata_list = []
for filepath in in_dir.glob("*.h5ad"):
    sample = filepath.stem.split("_")[2]
    print(f'Processing sample {sample}..')
    
    ad = sc.read_h5ad(filepath)
    
    # Filtering singlets only - OPTIONAL 
    # ad = ad[ad.obs["scDblFinder.class"] == "singlet"].copy()

    # Reading in file
    ad.obs['sample'] = sample
    
    ad.obs['treatment'] = np.where(ad.obs['sample'].isin(['BC3', 'BC9', 'BC15']), 'CORT', 'OIL')

    adata_list.append(ad)

    print(f'Added {filepath} to the list!')

# ========== CONCATENATING ANNDATA OBJECTS ==========
adata_all = anndata.concat(
    adata_list, 
    join='inner',  # Use only shared genes
    merge='same',  # Keep .obs/.var columns only if all match
    index_unique=None
)

print(adata_all.obs.head())

adata_all.obs_names_make_unique() # some barcodes are shared across samples 

adata_all.write(os.path.join(out_dir,f'Merged_n={adata_all.n_obs}_all_metadata_slide_tags.h5ad'))

