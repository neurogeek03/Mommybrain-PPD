from functions import filter_1to1_mouse_orthologs
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import numpy as np


stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# Custom paths
data_dir = project_path / 'data' 
slide_tags_merged = data_dir / 'add_mouse_orthologs' / 'MT_removed_merged_filtered_129493_mincells_10_in_2_samples_slide_tags.h5ad'
rat_to_mouse_list = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/rat_to_mouse_filtered.csv')

# =================== PARAMS ===================
sample = 'BC13'
mouse_id_column = 'mouse_ID'

# =================== INPUT ===================
adata_all = sc.read_h5ad(slide_tags_merged)

# TODO need to FIX OBS 
correct_list = pd.read_csv(rat_to_mouse_list)
print(correct_list.head())

# remove unnecessary .var columns
columns_to_remove = [
   'mouse_gene_stable_ID',
   'mouse_gene_name',
   'mouse_protein_transcript_stable_ID',
   'mouse_chromosome_name',
   'mouse_orthology_confidence_0_1',
   'rat_gene_name',
   'gene_symbols'
   ]

existing_columns_to_drop = [col for col in columns_to_remove if 
                            col in adata_all.var.columns]
adata_all.var.drop(columns=existing_columns_to_drop, inplace=True)
print(f"\nSuccessfully removed {len(existing_columns_to_drop)} columns.")

# collapse to mouse gene ids only 
collapsed_adata_all = filter_1to1_mouse_orthologs(adata_all, mouse_id_col=mouse_id_column)

print('Object before coolapsing:')
print(adata_all)
print('Object AFTER collapsing - var names are mouse gene ids:')
print(collapsed_adata_all)

# =================== NORMALIZATION ===================
print('checking if there are any negative values...')
contains_negative = collapsed_adata_all.X.min() < 0
print(contains_negative)
print('normalizing..')
collapsed_adata_all
sc.pp.normalize_total(collapsed_adata_all, target_sum=1e4)
print('=====================================')
print('checking if there are any negative values...')
contains_negative = collapsed_adata_all.X.min() < 0
print(contains_negative)
print(collapsed_adata_all.X[:5, :5].toarray())
# =================== SAVE LARGE OBJECT ===================
collapsed_path = data_dir / 'collapsed_MT_removed_merged_filtered_129493_mincells_10_in_2_samples_slide_tags.h5ad'
collapsed_adata_all.write(collapsed_path)

# =================== SUBSET TO 1 SAMPLE ===================
adata_subset = collapsed_adata_all[collapsed_adata_all.obs['sample'] == sample].copy()
print(adata_subset)

# =================== DATA VALIDATION ===================
print("\n--- Starting Data Validation ---")

# Filter out genes with no symbol before validation
adata_subset = adata_subset[:, ~adata_subset.var['gene_symbol'].isnull()].copy()
print(f"  - Removed genes with no symbol. New shape: {adata_subset.shape}")

# Replace the index (ENSEMBL IDs) with gene symbols
adata_subset.var.index = adata_subset.var['gene_symbol'].astype(str)
print("  - Replaced index with gene symbols.")

# Ensure gene symbols are unique, appending numbers if necessary
adata_subset.var_names_make_unique()
print("  - Ensured gene names are unique.")
print("  - Preview of new .var (first 5):")
print(adata_subset.var.head())
# Temporarily create DataFrames for validation
validation_expr_df = adata_subset.to_df()
validation_meta_df = adata_subset.obs

# =================== OUTPUT ===================
expression_df = adata_subset.to_df()
transposed_expression_df = expression_df.T

# adding the cell subclass label as the last column
expression_df['cell_subclass'] = adata_subset.obs['subclass_name']
metadata_df = adata_subset.obs

# # =================== Saving ===================
csv_path = output_base / f"pos_{sample}_expression_matrix_cell_subclass.csv"
# index=True saves cell IDs as the first column
expression_df.to_csv(csv_path, index=True) 
meta_path = output_base / f"pos_{sample}_metadata.csv"
metadata_df.to_csv(meta_path, index=True)
