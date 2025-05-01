"""
Title:        Cellbender Output to .h5ad file 
Description:  Creating an anndata object from the gene by cell matrix produced by cellbender 
Author:       Maria Eleni Fafouti 
Date:         30-04-2025
"""

# bash
# conda activate seurat_env

# ========== IMPORTS ==========
import h5py
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os

# ========== DEFINING ESSENTIAL PATHS ==========
data_path = "/project/s/shreejoy/hudsonhu/CB_RatDataOnlyRun"
project_path = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags"
output_path = os.path.join(project_path, "Post_bender")
ortholog_df_path = os.path.join(project_path, "Gene_lists", "rat_mouse_genes.txt")
email = "mariaeleni.fafouti@mail.utoronto.ca"

# ========== MAPPING RAT GENE IDS TO MOUSE ORTHOLOGS ==========

# ===== Step 1: Filter the gene list
ortholog_df = pd.read_csv(ortholog_df_path, sep=",")
initial_rows = len(ortholog_df)
print(f'The initial rows of the reference are: {initial_rows}')
# Fiiltering the list 
ortholog_highconf = ortholog_df[ortholog_df["mouse_orthology_confidence_0_1"] == 1.0]
after_conf_filter = len(ortholog_highconf)
print(f"{after_conf_filter} have 1 -1 orthology to mouse.")
# Removing duplicates by rat ID (keeping first occurrence) 
ortholog_unique = ortholog_highconf.drop_duplicates(subset="rat_ID")
after_dedup = len(ortholog_unique)
print(f"{after_dedup} rows remain after removing duplicate rat_IDs.")

# ========== DEFINING SAMPLE AND TREATMENT LISTS ==========
# sample_list = ["BC13"]
sample_list = ["BC13", "BC14", "BC15", "BC28", "BC3", "BC9"]
cort_samples = ['BC3', 'BC9', 'BC15']

# ===== Step 2: Read the gene list for each sample, left join to the list
# Goal: create gene metadata of the SAME length and order
for sample in sample_list:
    print(f'Sample {sample} is being processed')
    h5_file = os.path.join(data_path, f"{sample}", "output_file_filtered.h5")

    # Loading metadata from h5 file 
    with h5py.File(h5_file, "r") as f:
        # Genes
        gene_ids = f["matrix/features/id"][:].astype(str)
        gene_names = f["matrix/features/name"][:].astype(str)

        # Barcodes
        barcodes = f["matrix/barcodes"][:].astype(str)

        # Sparse matrix components
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = f["matrix/shape"][:]

    # Building a dataframe for var 
    var_df = pd.DataFrame({
        "ensembl_id": gene_ids,
        "gene_symbol": gene_names
    }, index=gene_ids)

    # Left join by rat_ID
    merged = var_df.merge(ortholog_unique, how="left", left_on="ensembl_id", right_on="rat_ID")

    # Optional preview:
    print(f'Sample {sample} has the following list for var metadata')
    pd.set_option('display.max_columns', None)  # Show all columns
    print(merged.head())
    print(merged.notna().sum())  # See how many genes didn't match

    # Saving how many columns are full in the list for each sample 
    # Ensure the directory exists before saving
    var_meta_dir = os.path.join(output_path, f"{sample}")  # Sample directory
    var_meta_path = os.path.join(var_meta_dir, "var_meta_columns_qc.txt")  # Full file path
    os.makedirs(var_meta_dir, exist_ok=True)  # Create the directory (not just the parent)

    # Now save the file
    merged.notna().sum().to_csv(var_meta_path, sep="\t")

    # ===== Step 3: Create the Anndata object
    # Reconstruct CSC matrix: genes (rows) × cells (columns)
    X = sp.csc_matrix((data, indices, indptr), shape=shape)

    # Create AnnData object with Ensembl IDs as var_names (genes)
    adata = ad.AnnData(X.transpose())  # Now: cells × genes
    adata.var_names = gene_ids # gene labels 
    adata.obs_names = barcodes # cell labels 

    # Defining the treatment values, which will be saved in obs
    treatment = "CORT" if {sample} in cort_samples else "OIL"

    # Create .obs with relevant metadata
    adata.obs = pd.DataFrame({
        'barcode': barcodes,
        'sample': f"{sample}",
        'treatment': treatment
    }, index=barcodes)

    # ===== Step 4: Update the .var dataframe in AnnData
    # Now assign the merged data (with orthologs and gene symbols) to .var in the AnnData object
    # Ensure that merged's index matches the gene_ids (same order as adata.var_names)
    merged = merged.set_index('ensembl_id')  # Ensure the merge index is properly aligned
    adata.var = merged  # Assign merged data to the .var slot

    # Optional: Store gene symbols as metadata
    adata.var["gene_symbols"] = gene_names

    # Save to file if needed
    adata.write_h5ad(os.path.join(output_path, f"{sample}", f"converted_ann_data_{sample}.h5ad"))



# #TESTING to preview one sample's gene names \

# h5_file = os.path.join(data_dir, "BC13", "output_file_filtered.h5")

# # Function to preview file - to identify which specific datasets you wish to have in your anndata obj
# def print_structure(name, obj):
#     if isinstance(obj, h5py.Dataset):
#         print(f"{name}: dataset, shape = {obj.shape}, dtype = {obj.dtype}")
#     elif isinstance(obj, h5py.Group):
#         print(f"{name}: group")

# # Preview structure 
# with h5py.File(h5_file, "r") as f:
#     f.visititems(print_structure)

# with h5py.File(h5_file, "r") as f:
#     print("Keys under metadata:", list(f["metadata"].keys()))
#     # Load and preview gene names
#     gene_names = f["metadata/gene_names"][:].astype(str)
#     print(f"Total gene names: {len(gene_names)}")
#     print("First 20 gene names:", gene_names[:20])
