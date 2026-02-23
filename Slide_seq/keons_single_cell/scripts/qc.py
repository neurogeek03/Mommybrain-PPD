import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import seaborn as sns

from functions import plot_cell_counts, classify_cell_type

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# data
data_path = project_path / 'data' / 'All_RCTD_types_singlet_score_0_slide_seq_15.h5ad'

# =================== PARAMS ===================
qc_png = 'qc_metrics.png'
out_png = output_base / qc_png
fig_dir = output_base / 'figures'
fig_dir.mkdir(exist_ok=True, parents=True)
sc.settings.figdir = fig_dir

obs_column = 'RCTD_first_type_rat'
# =================== FUNCTIONS ===================

# =================== INPUT ===================
adata = sc.read_h5ad(data_path)

# sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
# sc.pl.violin(
#     adata, 
#     ['n_genes_by_counts', 'total_counts'], 
#     multi_panel=True, 
#     save="_qc_metrics.png"  # The prefix 'violin' is added automatically
# )

# =================== ADDING METADATA COLUMNS ===================

# column for first type NN
adata.obs['is_NN'] = adata.obs[obs_column].str.endswith('NN', na=False)
print(adata.obs['is_NN'].value_counts())
print(adata.obs[[obs_column, 'is_NN']].head(10))

# =================== SPLITTING DATA ===================
adata_non_neurons = adata[adata.obs['is_NN'] == True].copy()
adata_neurons = adata[adata.obs['is_NN'] == False].copy()

print(f"Original size: {adata.n_obs} beads")
print(f"Non-Neuron object: {adata_non_neurons.n_obs} beads")
print(f"Neuron object: {adata_neurons.n_obs} beads")

# adata_nn_subset = adata_non_neurons[adata_non_neurons.obs['RCTD_spot_class_rat'] =='singlet'].copy()
# adata_nn_subset = adata_nn_subset[adata_nn_subset.obs['RCTD_singlet_score_rat'] > 330].copy()

# adata_nn_subset_path = output_base / f'adata_nn_subset_{adata_nn_subset.n_obs}_singlets_score_330.h5ad'
# adata_nn_subset.write(adata_nn_subset_path)

########
adata_neuron_subset = adata_neurons[adata_neurons.obs['RCTD_singlet_score_rat'] > 330].copy()

adata_neuron_subset_path = output_base/f'adata_nn_subset_{adata_neuron_subset.n_obs}_singlets_score_330.h5ad'
adata_neuron_subset.write(adata_neuron_subset_path)

# sc.pp.calculate_qc_metrics(adata_neurons, percent_top=None, log1p=False, inplace=True)
# sc.pl.violin(
#     adata_neurons, 
#     ['n_genes_by_counts', 'total_counts'], 
#     multi_panel=True, 
#     save="neurons_qc_metrics.png"  # The prefix 'violin' is added automatically
# )

# sc.pp.calculate_qc_metrics(adata_non_neurons, percent_top=None, log1p=False, inplace=True)
# sc.pl.violin(
#     adata_non_neurons, 
#     ['n_genes_by_counts', 'total_counts'], 
#     multi_panel=True, 
#     save="nn_qc_metrics.png"  # The prefix 'violin' is added automatically
# )

# # =================== COUNTS PLOTS ===================
# # 1. Plot for Non-Neurons (The "NN" group)
# plot_cell_counts(adata_non_neurons,obs_column,"Cell Type Counts: Non-Neurons (NN)", "counts_non_neurons.png", use_log=True)

# # 2. Plot for Neurons (The non-NN group)
# plot_cell_counts(adata_neurons,obs_column, "Cell Type Counts: Neurons", "counts_neurons.png", use_log=True)

# =================== OUTPUT ===================

# ========== SCRATCH

# # Renaming the var part
# adata.var_names = adata.var_names.str.split("-").str[-1]
# print(adata.var_names[:10])
# is_unique = adata.var_names.is_unique
# print(is_unique) # excellent
# print(adata.var_names[:10])

# # checking for mt genes in var 
# cols_to_check = ['name', 'gene_symbol', 'gene_id']

# for col in cols_to_check:
#     # Check if the column exists first
#     if col in adata.var.columns:
#         # Regex '^(([Mm][Tt])-)' looks for MT-, mt-, or Mt- at the start of the string
#         is_mito = adata.var[col].str.contains('^mt-', case=False, na=False)
#         count = is_mito.sum()
#         print(f"Column '{col}': found {count} mitochondrial genes.")
        
#         # Optional: Print the first few matches to verify
#         if count > 0:
#             print(f"  Examples: {adata.var[col][is_mito].head(3).values}")