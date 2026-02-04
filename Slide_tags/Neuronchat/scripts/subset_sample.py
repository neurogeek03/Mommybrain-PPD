from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import numpy as np
import argparse 
from functions import find_file_by_parts

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

#TODO edit input file names to be dependent on the sample arg in the R script 

# =================== ARGS ===================
parser = argparse.ArgumentParser(description="Preprocessing input for neuronchat.")
parser.add_argument("-i", "--input", help="Path to your input file with all samples" )
parser.add_argument("-s", "--sample", help="The sample you wish to run" )
parser.add_argument("-o", "--output_base", help="Folder where outputs are stored " )

args = parser.parse_args()

# =================== PARAMS ===================
sample = args.sample 


# =================== PATHS ===================
# project_path = Path.cwd().parents[0]
# print(f"current working directory: {project_path}")
input_dir = Path(args.input)
output_base = Path(args.output_base)
output_base.mkdir(exist_ok=True, parents=True)

#FIXME change this relative to argparse
files_matching = find_file_by_parts(input_dir, "collapsed", ".h5ad")
adata_path = files_matching[0]
# merged_samples_path = input_dir / 'collapsed_MT_removed_merged_filtered_129493_mincells_10_in_2_samples_slide_tags.h5ad'

# =================== INPUT ===================
collapsed_adata_all = sc.read_h5ad(adata_path)

# =================== SUBSET TO 1 SAMPLE ===================
adata_subset = collapsed_adata_all[collapsed_adata_all.obs['sample'] == sample].copy()
print(adata_subset)

# # =================== DATA VALIDATION ===================
# print("\n--- Starting Data Validation ---")

# # Filter out genes with no symbol before validation
# adata_subset = adata_subset[:, ~adata_subset.var['gene_symbol'].isnull()].copy()
# print(f"  - Removed genes with no symbol. New shape: {adata_subset.shape}")

# # Replace the index (ENSEMBL IDs) with gene symbols
# adata_subset.var.index = adata_subset.var['gene_symbol'].astype(str)
# print("  - Replaced index with gene symbols.")

# # Ensure gene symbols are unique, appending numbers if necessary
# adata_subset.var_names_make_unique()
# print("  - Ensured gene names are unique.")
# print("  - Preview of new .var (first 5):")
# print(adata_subset.var.head())
# # Temporarily create DataFrames for validation
# validation_expr_df = adata_subset.to_df()
# validation_meta_df = adata_subset.obs

# =================== OUTPUT ===================
expression_df = adata_subset.to_df()
transposed_expression_df = expression_df.T

# adding the cell subclass label as the last column
#FIXME here we can also have class labels. So you have to pick which one you want every time. 
expression_df['cell_subclass'] = adata_subset.obs['subclass_name']
metadata_df = adata_subset.obs

# # =================== Saving ===================
csv_path = output_base / f"{sample}_expression_matrix_cell_subclass.csv"
# index=True saves cell IDs as the first column
expression_df.to_csv(csv_path, index=True) 
meta_path = output_base / f"{sample}_metadata.csv"
metadata_df.to_csv(meta_path, index=True)