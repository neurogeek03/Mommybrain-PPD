import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import harmonypy as hm
import numpy as np
from pathlib import Path

input_dir = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Filtering')
adata_dir = input_dir / 'NEW_list_merged_filtered' 
adata_umap = adata_dir /'umap_filtered_150725_slide_tags.h5ad'

adata = sc.read_h5ad(adata_umap)

print('done reading data!')

# ========== RUNNING HARMONY ==========
sce.pp.harmony_integrate(adata, key="sample", max_iter_harmony=20)  # Creates `X_pca_harmony` in adata.obsm

# ========== USING HARMONY-ADJUSTED PCs DOWNSTREAM: NEIGHBORS & CLUSTERING ==========
print('finiding neighbors...')
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
print('running UMAP...')
sc.tl.umap(adata)
print('running leiden clustering...')
sc.tl.leiden(adata)

out_path = adata_dir / f"harmony_umap_filtered_{adata.n_obs}_slide_tags.h5ad"
adata.write(out_path)