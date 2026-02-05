from functions import find_file_by_parts, filter_1to1_mouse_orthologs, collapse_by_gene_symbol
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import numpy as np
import argparse 

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# Custom paths
data_dir = project_path / 'data' / 'raw'
files_matching = find_file_by_parts(data_dir, "DE", ".h5ad")
adata_path = files_matching[0]

# slide_tags_merged = data_dir / 'MT_removed_merged_filtered_129493_mincells_10_in_2_samples_slide_tags.h5ad'
rat_to_mouse_list = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/rat_to_mouse_filtered.csv')

# =================== PARAMS ===================
mouse_id_column = 'Mouse gene stable ID'
mouse_symbol_column = 'Mouse gene name'

# =================== INPUT ===================
adata_all = sc.read_h5ad(adata_path)

correct_list = pd.read_csv(rat_to_mouse_list)
print(correct_list.head())

# =================== ADD MOUSE METADATA ===================
adata_all.var = pd.merge(
    adata_all.var,
    correct_list,
    left_index=True,
    right_on='Gene stable ID',
    how='left'
)
adata_all.var.set_index(adata_all.var_names, inplace=True)

print(adata_all)

# collapse to mouse gene ids only 
collapsed_adata_all = filter_1to1_mouse_orthologs(adata_all, mouse_id_col=mouse_id_column)
print('Object before coolapsing:')
print(adata_all)
print('Object AFTER collapsing - var names are mouse gene ids:')
print(collapsed_adata_all)


# Collapse to gene names 
collapsed_by_symbol = collapse_by_gene_symbol(collapsed_adata_all, gene_symbol_col=mouse_symbol_column)
print('Object before coolapsing BY SYMBOL:')
print(collapsed_adata_all)
print('Object AFTER collapsing - var names are mouse gene SYMBOLS:')
print(collapsed_by_symbol)


# =================== SAVE RAW COUNTS AS LAYER ===================
collapsed_by_symbol.layers["counts"] = collapsed_by_symbol.X.copy()

# =================== NORMALIZATION ===================
print('checking if there are any negative values...')
contains_negative = collapsed_by_symbol.X.min() < 0
print(contains_negative)
print('normalizing..')
collapsed_by_symbol
sc.pp.normalize_total(collapsed_by_symbol, target_sum=1e4)
print('=====================================')
print('checking if there are any negative values...')
contains_negative = collapsed_by_symbol.X.min() < 0
print(contains_negative)
# print(collapsed_by_symbol.X[:5, :5].toarray())
print(collapsed_by_symbol.X[:5, :5])

# =================== SAVE LARGE OBJECT ===================
collapsed_path = data_dir / 'layers_collapsed_DE_symbols.h5ad'
collapsed_by_symbol.write(collapsed_path)


# SCRATCH
# # remove unnecessary .var columns
# columns_to_remove = [
#    'mouse_gene_stable_ID',
#    'mouse_gene_name',
#    'mouse_protein_transcript_stable_ID',
#    'mouse_chromosome_name',
#    'mouse_orthology_confidence_0_1',
#    'rat_gene_name',
#    'gene_symbols'
#    ]

# existing_columns_to_drop = [col for col in columns_to_remove if 
#                             col in adata_all.var.columns]
# adata_all.var.drop(columns=existing_columns_to_drop, inplace=True)
# print(f"\nSuccessfully removed {len(existing_columns_to_drop)} columns.")