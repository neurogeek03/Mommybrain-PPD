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
subset_dir = output_base / 'subset_objects'
subset_dir.mkdir(exist_ok=True, parents=True)

# data
data_path = project_path / 'data' / 'All_RCTD_types_singlet_score_0_slide_seq_15.h5ad'

# =================== PARAMS ===================
fig_dir = output_base / 'figures'
fig_dir.mkdir(exist_ok=True, parents=True)
sc.settings.figdir = fig_dir

obs_column = 'RCTD_first_type_rat'

# =================== INPUT ===================
adata = sc.read_h5ad(data_path)

# Keep only the prefix before '-' in the var 'name' column
adata.var['name'] = adata.var['name'].str.split('-').str[0]

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

# =================== FILTER & SAVE SUBSETS ===================
adata_nn_subset = adata_non_neurons[adata_non_neurons.obs['RCTD_spot_class_rat'] == 'singlet'].copy()
adata_nn_subset = adata_nn_subset[adata_nn_subset.obs['RCTD_singlet_score_rat'] > 330].copy()
print(f"Non-Neuron after filtering: {adata_nn_subset.n_obs} beads")

adata_neuron_subset = adata_neurons[adata_neurons.obs['RCTD_singlet_score_rat'] > 330].copy()
print(f"Neuron after filtering: {adata_neuron_subset.n_obs} beads")

adata_nn_subset.write(subset_dir / f'adata_non_neurons_{adata_nn_subset.n_obs}.h5ad')
adata_neuron_subset.write(subset_dir / f'adata_neurons_{adata_neuron_subset.n_obs}.h5ad')
print(f"Subsets saved to {subset_dir}")

# =================== NON-NEURON PROCESSING ===================

nn_subclasses = list(adata_non_neurons.obs['RCTD_first_type_rat'].cat.categories)
print(nn_subclasses)

# Subclass → Class mapping (Allen Brain CCN20230722 taxonomy)
# Source: CLASS_exlcusion_baseline_csv_mapping_output.csv
nn_subclass_to_class = {
    '318_Astro_NT_NN':       '30 Astro-Epen',
    '319_Astro_TE_NN':       '30 Astro-Epen',
    '320_Astro_OLF_NN':      '30 Astro-Epen',
    '321_Astroependymal_NN': '30 Astro-Epen',
    '322_Tanycyte_NN':       '30 Astro-Epen',
    '323_Ependymal_NN':      '30 Astro-Epen',
    '325_CHOR_NN':           '30 Astro-Epen',
    '326_OPC_NN':            '31 OPC-Oligo',
    '327_Oligo_NN':          '31 OPC-Oligo',
    '328_OEC_NN':            '32 OEC',
    '329_ABC_NN':            '33 Vascular',
    '330_VLMC_NN':           '33 Vascular',
    '331_Peri_NN':           '33 Vascular',
    '332_SMC_NN':            '33 Vascular',
    '333_Endo_NN':           '33 Vascular',
    '334_Microglia_NN':      '34 Immune',
    '335_BAM_NN':            '34 Immune',
    '337_DC_NN':             '34 Immune',
    '338_Lymphoid_NN':       '34 Immune',
}

adata_non_neurons.obs['allen_class'] = adata_non_neurons.obs[obs_column].map(nn_subclass_to_class)
print(adata_non_neurons.obs['allen_class'].value_counts())
print(adata_non_neurons.obs['allen_class'].isna().sum(), "unmapped beads")

adata_nn_subset = adata_non_neurons[adata_non_neurons.obs['RCTD_spot_class_rat'] == 'singlet'].copy()

adata_nn_subset_path = output_base / f'class_singlets_broader_adata_nn_subset_{adata_nn_subset.n_obs}_score_330.h5ad'
adata_nn_subset.write(adata_nn_subset_path)