import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import warnings
from pathlib import Path
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod
from sccoda.util import data_visualization as viz

# params 
comparison_folder = 'caudal_cort_vs_non_cort' # change name of folder at each iteration
pd.set_option('display.max_columns', None)

# ============== paths ==============
project_path = Path.cwd()
data_dir = project_path / 'data' 
data_path = data_dir / 'cell_counts_per_type_1M_slide_seq.csv'
cov_path = data_dir / 'covariates_caudal.csv'
out_dir = project_path / 'out' / 'slide_seq' / comparison_folder
out_dir.mkdir(exist_ok=True, parents=True)

# # ============== Load data from scanpy ==============
# adata = ad.read_h5ad(data_path)
# create covariate df 
# cov_df = adata.obs[['sample','treatment']].copy().drop_duplicates().set_index('sample')
# print(cov_df.head())
# data_scanpy_1 = dat.from_scanpy(
#     adata,
#     cell_type_identifier="subclass_name",
#     sample_identifier="sample",
#     covariate_df=cov_df
# )
# print(data_scanpy_1.to_df().head())
# print(data_scanpy_1)
# cell_counts_obj = data_scanpy_1

# ============== Load data from manually-made csv file ==============
# read in custom covariate df
cov_df = pd.read_csv(cov_path)
print(cov_df.head())
print('Samples included in the covariates df:')
print(cov_df['sample'])

# read in cell counts for all data 
cell_counts = pd.read_csv(data_path)

# subset cell counts based on covariate index
cell_counts = cell_counts[cell_counts["sample"].isin(cov_df["sample"])]
print('Samples included in the cell counts df:')
print(cell_counts['sample'])

# Filtering based on neuron type
cell_counts=cell_counts.loc[:, ~cell_counts.columns.str.contains('NN$')]
#print(cell_counts.head())

# use pandas df to create sccoda anndata obj
data_slide_seq = dat.from_pandas(cell_counts, covariate_columns=["sample"])
print(data_slide_seq)
print(data_slide_seq.to_df().head())
data_slide_seq.obs = data_slide_seq.obs.merge(
    cov_df,
    on="sample",
    how="left"  # keep all cells in adata.obs
)
print(data_slide_seq.X)
print(data_slide_seq.obs)

cell_counts_obj = data_slide_seq

# ax = viz.boxplots(cell_counts_obj, feature_name="condition")
# fig = ax.get_figure()
# fig.set_size_inches(25, 6)  # for example, 10x6 inches
# ax.tick_params(axis='x', labelsize=8, rotation=90) 
# fig.tight_layout()
# fig.savefig(out_dir / "slide_seq_boxplot.png", dpi=300, bbox_inches="tight")

# RUNNING
model = mod.CompositionalAnalysis(cell_counts_obj, formula = "condition")
results = model.sample_hmc()  # runs the Bayesian inference
print(results.summary())
print(results.effect_df.columns)

results.effect_df.to_csv(out_dir / "sccoda_effects.csv")
results.intercept_df.to_csv(out_dir / "sccoda_intercepts.csv")

