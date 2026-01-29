from functions import filter_1to1_mouse_orthologs
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import numpy as np
import argparse 

#TODO edit input file names to be dependent on the sample arg in the R script 

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
mouse_id_column = 'mouse_ID'

# =================== INPUT ===================
adata_all = sc.read_h5ad(slide_tags_merged)

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


