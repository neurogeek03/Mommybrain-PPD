import scanpy as sc
import pandas as pd

# Load the AnnData object
adata = sc.read_h5ad("/scratch/mfafouti/Mommybrain/Slide_seq/Integration/umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad")  # replace with your file path

# Make sure 'sample' and 'RCTD_first_type_mouse' are in adata.obs
# Count cells per sample and RCTD_first_type_mouse
counts = adata.obs.groupby(['sample', 'RCTD_first_type_mouse']).size().reset_index(name='count')

# Keep only the RCTD_first_type_mouse where each sample has at least 10 cells
# For that, sum counts per sample per type and filter
counts_filtered = counts[counts['count'] >= 10]

# Optional: if you want to keep all RCTD_first_type_mouse even if count <10 but only show types with >=10 somewhere else
# counts_filtered = counts.copy() # for the table you showed

# Save as CSV
counts_filtered.to_csv("RCTD_counts_filtered.csv", index=False)
