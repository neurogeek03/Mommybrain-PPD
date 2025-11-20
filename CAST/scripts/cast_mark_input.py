import os
from functions import preview_adata_X 
import numpy as np
import anndata as ad
import scanpy as sc
# import argparse
import json
from pathlib import Path
import warnings

print('imports done')

# Argumet parsing
# parser = argparse.ArgumentParser()
# parser.add_argument("-sample", type=str, required=True)
# args = parser.parse_args()

# config 
# sample = args.sample
warnings.filterwarnings("ignore")

# paths 
work_dir = Path.cwd()
data_path = work_dir / 'data' 
output_path = work_dir / 'out'
os.makedirs(output_path, exist_ok=True)

# per sample 
merged_adata_path = data_path / 'RAW_slideseq_1054147_same_as_umap.h5ad'
adata = ad.read_h5ad(merged_adata_path)
print('sample read')

preview = preview_adata_X(adata, 20)
print(preview)

# normalize raw gene expression data
adata.layers['norm_1e4'] = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)['X'].toarray() # we use normalized counts for each cell as input gene expression

# Extract information per sample 
samples = np.unique(adata.obs['sample']) # used samples in adata
coords_raw = {sample_t: adata.obsm["X_spatial"][adata.obs['sample'] == sample_t, :]
    for sample_t in samples}
exp_dict = {sample_t: adata[adata.obs['sample'] == sample_t].layers['norm_1e4'] 
            for sample_t in samples}

with open(output_path, "w") as f:
    json.dump({"coords_dict": coords_raw, "exp_dict": exp_dict}, f)
