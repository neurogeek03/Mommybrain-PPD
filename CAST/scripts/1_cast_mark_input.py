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
adata_dir = work_dir / 'out' / 'objects'

# ============================= PROCESS MERGED =============================
merged_adata_path = adata_dir /'B01_B42_subset.h5ad'
merged = ad.read_h5ad(merged_adata_path)
print('OBJECT read')

# normalize raw gene expression data
merged.layers['norm_1e4'] = sc.pp.normalize_total(merged, target_sum=1e4, inplace=False)['X'].toarray() # we use normalized counts for each cell as input gene expression
print('data norm done')
# Extract information per sample 
samples = np.unique(merged.obs['sample']) # used samples in adata

# Modified code to avoid memory overload
coords_dict = {
    s: merged.obsm["X_spatial"][merged.obs['sample'] == s, :]
    for s in samples
} # WORKS UP TO HERE!!!
exp_dict = {
    s: merged.layers["norm_1e4"][merged.obs['sample'] == s, :]
    for s in samples
}

print('data dict creation done')


coords_clean = {k: v.tolist() for k, v in coords_dict.items()}
exp_clean = {k: v.tolist() for k, v in exp_dict.items()}
print('data dict cleaning done')
with open(output_path, "w") as f:
    json.dump({"coords_dict": coords_clean, "exp_dict": exp_clean}, f)
print(f'data saved at {output_path}')