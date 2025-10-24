# export NUMBA_CACHE_DIR=/scratch/mfafouti/Mommybrain/Slide_tags/Proportion_analysis/numba_cache

import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import warnings
from pathlib import Path

from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

project_path = Path.cwd().parents[0]
data_dir = project_path / 'data' 
data_path = data_dir / 'PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'

adata = ad.read_h5ad(data_path)

# create covariate df 
cov_df = adata.obs[['sample','treatment']].copy()
print(cov_df.head())

# Load data
data_scanpy_1 = dat.from_scanpy(
    adata,
    cell_type_identifier="subclass_name",
    sample_identifier="sample",
    covariate_df=cov_df
)
print(data_scanpy_1)