from functions import filter_1to1_mouse_orthologs
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc


stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# Custom paths
data_dir = project_path / 'data' 
slide_tags_merged = data_dir / 'PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'
rat_to_mouse_list = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/rat_to_mouse_filtered.csv')

# =================== PARAMS ===================
sample = 'BC13'
mouse_id_column = 'mouse_gene_stable_ID'

# =================== INPUT ===================
adata_all = sc.read_h5ad(slide_tags_merged)

# TODO need to FIX OBS 


# collapse to mouse gene ids only 
collapsed_adata_all = filter_1to1_mouse_orthologs(adata_all, mouse_id_col=mouse_id_column)

print('Object before coolapsing:')
print(adata_all)
print('Object AFTER collapsing - var names are mouse gene ids:')
print(collapsed_adata_all)

# collapsed_path = data_dir / 'collapsed_PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'
# collapsed_adata_all.write(collapsed_path)

# =================== SUBSET TO 1 SAMPLE ===================
adata_subset = adata_all[adata_all.obs['sample'] == sample].copy()
print(adata_subset)

# =================== OUTPUT ===================
expression_df = adata_subset.to_df()
# index=True saves cell IDs as the first column
csv_path = output_base / "expression_matrix.csv"
expression_df.to_csv(csv_path, index=True) 

metadata_df = adata_subset.obs
meta_path = output_base /"metadata.csv"
metadata_df.to_csv(meta_path, index=True)
