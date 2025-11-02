# export NUMBA_CACHE_DIR=/scratch/mfafouti/Mommybrain/Slide_tags/Proportion_analysis/numba_cache

import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import warnings
from pathlib import Path

from sccoda.util import cell_composition_data as dat
# from sccoda.model import CompositionalAnalysis
from sccoda.util import comp_ana as mod
from sccoda.util import data_visualization as viz

project_path = Path.cwd()
data_dir = project_path / 'data' 
data_path = data_dir / 'PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'
out_dir = project_path / 'out'

# functions
from sccoda.util import data_visualization as viz
import matplotlib.pyplot as plt


adata = ad.read_h5ad(data_path)

# create covariate df 
cov_df = adata.obs[['sample','treatment']].copy().drop_duplicates().set_index('sample')
print(cov_df.head())

# Load data
data_scanpy_1 = dat.from_scanpy(
    adata,
    cell_type_identifier="subclass_name",
    sample_identifier="sample",
    covariate_df=cov_df
)
print(data_scanpy_1)

ax = viz.boxplots(data_scanpy_1, feature_name="treatment")
fig = ax.get_figure()
fig.set_size_inches(20, 6)  # for example, 10x6 inches
ax.tick_params(axis='x', labelsize=8, rotation=90) 
fig.tight_layout()
fig.savefig(out_dir / "test1_boxplot.png", dpi=300, bbox_inches="tight")

# RUNNING
model = mod.CompositionalAnalysis(data_scanpy_1, formula = "treatment")
results = model.sample_hmc()  # runs the Bayesian inference
print(results.summary())
print(results.effect_df.columns)

results.effect_df.to_csv(out_dir / "sccoda_effects.csv")
results.intercept_df.to_csv(out_dir / "sccoda_intercepts.csv")
