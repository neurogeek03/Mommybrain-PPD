"""
Title: Extract edgeR input from combined anndata object
Description: Pseudobulk for RCTD "singlets" data by summing raw counts for all samples in each cell type from a single AnnData file
Author:   Maria Eleni Fafouti 
Date: 17-07-2025
"""
import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from collections import defaultdict
import re

# === Input & output paths ===
current_path = '/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR'
adata_path = '/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects'
input_file = os.path.join(adata_path,  "All_RCTD_types_singlet_score_0_slide_seq_15.h5ad")
output_dir = os.path.join(current_path,"out","pseudobulk_outputs")
os.makedirs(output_dir, exist_ok=True)

celltype_col = "RCTD_first_type_rat"  
# === Functions ===
def clean_symbol(name):
    parts = name.split("-")
    if parts[0] == ".":   # case where it's just a dot
        return parts[-1]  # keep after the dash
    else:
        return parts[0]

# === Load the combined AnnData object ===
adata = sc.read_h5ad(input_file)

# =========== FIXING THE VAR SECTION ===========
# Remove the extra gene_id column (keep only the index)
adata.var = adata.var.drop(columns=["gene_id"])
# cleaning the gene_symbol column 
adata.var["gene_symbol"] = adata.var["name"].apply(clean_symbol).astype(str)
adata.write(os.path.join(adata_path,f'cleaned_adata_filtered_220626_10_in_any_2_samples_singlet_score_300.h5ad'))

# === Collect all genes ===
# Make gene_symbol the index / var_names as a string Index
adata.var_names = pd.Index(adata.var["gene_symbol"].values, dtype=str)

# Now enforce uniqueness safely
adata.var_names_make_unique()# ensures no duplicates
all_genes = adata.var_names.tolist()

# === Dictionary to hold counts for each cell type ===
pseudobulk_by_celltype = defaultdict(list)  # key = cell type, value = list of (sample_id, summed counts)

# === For each cell type, group by sample and sum counts ===
# iterate over cell types present in the dataset
celltypes = pd.unique(adata.obs[celltype_col].values)

for celltype in celltypes:
    # boolean mask for cells of this celltype
    mask_ct = adata.obs[celltype_col] == celltype
    # which samples actually contain this cell type
    samples_with_ct = adata.obs.loc[mask_ct, "sample"].unique()
    print(f"Processing celltype: {celltype}  (found in {len(samples_with_ct)} samples)")

    for samp in samples_with_ct:
        mask = mask_ct & (adata.obs["sample"] == samp) # for this specific cell type and sample
        X = adata[mask].X
        # robustly sum to a 1D numpy array whether sparse or dense:
        summed = X.sum(axis=0)
        # handle sparse/dense difference return types
        summed = np.asarray(summed).ravel()
        if summed.size != len(all_genes):
            raise RuntimeError(f"Length mismatch for {celltype}/{samp}: {summed.size} vs {len(all_genes)}")
        pseudobulk_by_celltype[celltype].append((samp, summed))
        print(f"  Added sample: {samp}")

print("Finished creating pseudobulk lists.")


# === Write one file per cell type ===
for celltype, sample_sums in pseudobulk_by_celltype.items():
    sample_ids, count_arrays = zip(*sample_sums)
    # create genes x samples matrix
    mat = np.column_stack(count_arrays)   # shape: (n_genes, n_samples_for_this_ct)
    df = pd.DataFrame(mat, index=all_genes, columns=sample_ids)
    # optionally cast to int if counts are integer-like:
    # df = df.astype(int)
    safe_name = str(celltype).replace("/", "_").replace(" ", "_")
    df.to_csv(f"{output_dir}/{safe_name}_counts.csv", index = True)