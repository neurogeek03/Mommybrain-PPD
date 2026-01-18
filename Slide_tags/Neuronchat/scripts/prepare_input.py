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
slide_tags_merged = data_dir / 'add_mouse_orthologs' / 'PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'
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

collapsed_path = data_dir / 'collapsed_PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'
collapsed_adata_all.write(collapsed_path)

# =================== SUBSET TO 1 SAMPLE ===================
adata_subset = collapsed_adata_all[collapsed_adata_all.obs['sample'] == sample].copy()
print(adata_subset)

# =================== DATA VALIDATION ===================
print("\n--- Starting Data Validation ---")
# Temporarily create DataFrames for validation
validation_expr_df = adata_subset.to_df()
validation_meta_df = adata_subset.obs

# 1. Verify metadata columns ('df_group' equivalent)
print("\n[1] Verifying metadata columns...")
required_meta_cols = ['class_label', 'subclass_label']
for col in required_meta_cols:
    if col in validation_meta_df.columns:
        print(f"  - Column '{col}' found. Unique values: {validation_meta_df[col].unique().tolist()}")
    else:
        print(f"  - WARNING: Column '{col}' NOT FOUND in metadata.")

# 2. Check for NA/NaN values
print("\n[2] Checking for NA/NaN values...")
if validation_expr_df.isnull().values.any():
    print("  - WARNING: NA/NaN values found in the expression data.")
else:
    print("  - OK: No NA/NaN values found in the expression data.")

if validation_meta_df.isnull().values.any():
    print("  - WARNING: NA/NaN values found in the metadata.")
else:
    print("  - OK: No NA/NaN values found in the metadata.")

# 3. Check for all-zero rows/columns in expression data
print("\n[3] Checking for all-zero rows/columns in expression data...")
# Note: This check is on the dense DataFrame, which can be memory intensive for large data
all_zero_rows = validation_expr_df.eq(0).all(axis=1).sum()
all_zero_cols = validation_expr_df.eq(0).all(axis=0).sum()
print(f"  - Found {all_zero_rows} rows (cells) that are all zeros.")
print(f"  - Found {all_zero_cols} columns (genes) that are all zeros.")
print("--- Data Validation Complete ---\n")

# =================== OUTPUT ===================
expression_df = adata_subset.to_df()
# adding the cell subclass label as the last column
expression_df['cell_subclass'] = adata_subset.obs['subclass_name']
metadata_df = adata_subset.obs

# =================== Saving ===================
csv_path = output_base / f"{sample}_expression_matrix_cell_subclass.csv"
# index=True saves cell IDs as the first column
expression_df.to_csv(csv_path, index=True) 
meta_path = output_base / f"{sample}_metadata.csv"
metadata_df.to_csv(meta_path, index=True)
