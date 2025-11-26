import os
from functions import preview_adata_X 
import numpy as np
import anndata as ad
import scanpy as sc
# import argparse
import json
from pathlib import Path
import gc

# paths 
work_dir = Path.cwd().parents[0]
print(f'Current working directroy {work_dir}')
data_path = work_dir / 'data'
query_path = data_path / 'query' / 'rat' 
# adata_dir = work_dir / 'out' / 'objects'

adata_dir = Path('/scratch/mfafouti/BANKSY/Banksy_py/data/query')

# ============================= LOAD COMBINED OBJECT AND SUBSET =============================
all_slideseq = query_path /'RAW_slideseq_1054147_same_as_umap.h5ad'
adata = ad.read_h5ad(all_slideseq)
print('OBJECT read')

samples_of_interest = ["B01"]
mask = adata.obs["sample"].isin(samples_of_interest)
adata_sub = adata[mask].copy()

print(adata_sub)

mask = adata_sub.obs["RCTD_spot_class_rat"] != "reject"
adata = adata_sub[mask].copy()

saving_path = adata_dir / 'B01_subset_from_RCTD_out.h5ad'
adata.write(saving_path)

# # ============================= MERGE 2 SAMPLES =============================
# # os.makedirs(output_path, exist_ok=True)
# # Loading & merging adata files 
# B01_adata = sc.read_h5ad(query_path / 'B01_with_RCTD_mouse.h5ad')
# B42_adata = sc.read_h5ad(query_path / 'B42_with_RCTD_mouse.h5ad')
# # Add sample identifiers
# B01_adata.obs['sample_id'] = 'B01'
# B42_adata.obs['sample_id'] = 'B42'
# print('data loading done')
# # checking if raw 
# # preview = preview_adata_X(B42_adata, 20) # yes 
# # Merge
# merged = ad.concat(
#     [B01_adata, B42_adata],
#     join='inner',
#     merge='same',
#     index_unique=None
# )
# print(merged)
# print('data merging done')
# adata_dir.mkdir(exist_ok=True, parents=True)
# adata_path = adata_dir / 'merged_B01_B42_raw.h5ad'
# merged.write(adata_path)

# del B01_adata, B42_adata
# gc.collect()
# print('data deleting done')