import os
from functions import preview_adata_X 
import numpy as np
import anndata as ad
import scanpy as sc
# import argparse
import json
from pathlib import Path
import warnings
import gc

print('imports done')

# Argumet parsing
# parser = argparse.ArgumentParser()
# parser.add_argument("-sample", type=str, required=True)
# args = parser.parse_args()

# config 
# sample = args.sample
warnings.filterwarnings("ignore")

# paths 
work_dir = Path.cwd().parents[0]
print(f'Current working directroy {work_dir}')
data_path = work_dir / 'data'
query_path = data_path / 'query' / 'rat' 
output_path = work_dir / 'out' / 'cast_mark_input.json'
# os.makedirs(output_path, exist_ok=True)

# Loading & merging adata files 
B01_adata = sc.read_h5ad(query_path / 'B01_with_RCTD_mouse.h5ad')
B42_adata = sc.read_h5ad(query_path / 'B42_with_RCTD_mouse.h5ad')

# Add sample identifiers
B01_adata.obs['sample_id'] = 'B01'
B42_adata.obs['sample_id'] = 'B42'

print('data loading done')
# checking if raw 
# preview = preview_adata_X(B42_adata, 20) # yes 

# Merge
merged = ad.concat(
    [B01_adata, B42_adata],
    join='inner',
    merge='same',
    index_unique=None
)

print(merged)
print('data merging done')
del B01_adata, B42_adata
gc.collect()
print('data deleting done')
# # per sample 
# # merged_adata_path = data_path / 'RAW_slideseq_1054147_same_as_umap.h5ad'
# # adata = ad.read_h5ad(merged_adata_path)
# print('sample read')

# normalize raw gene expression data
merged.layers['norm_1e4'] = sc.pp.normalize_total(merged, target_sum=1e4, inplace=False)['X'].toarray() # we use normalized counts for each cell as input gene expression
print('data norm done')
# Extract information per sample 
samples = np.unique(merged.obs['sample_id']) # used samples in adata
coords_raw = {sample_t: merged.obsm["X_spatial"][merged.obs['sample_id'] == sample_t, :] for sample_t in samples}
exp_dict = {sample_t: merged[merged.obs['sample_id'] == sample_t].layers['norm_1e4'] for sample_t in samples}
print('data dict creation done')
coords_clean = {k: v.tolist() for k, v in coords_raw.items()}
exp_clean = {k: v.tolist() for k, v in exp_dict.items()}
print('data dict cleaning done')
with open(output_path, "w") as f:
    json.dump({"coords_dict": coords_clean, "exp_dict": exp_clean}, f)
print(f'data saved at {output_path}')