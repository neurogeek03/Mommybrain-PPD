import anndata as ad
from pathlib import Path

wd = Path.cwd().parents[0]
print(f'Current project directory {wd}')
object_dir = wd / 'FINAL_run_newgenelist' / 'objects'

# Suppose you have two AnnData objects
adata = ad.read_h5ad("/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects/All_RCTD_types_singlet_score_0_slide_seq_15.h5ad")
adata_ref = ad.read_h5ad("/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects/1054147_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad")

# Extract the barcodes (cell IDs) from the reference
ref_barcodes = set(adata_ref.obs_names)

# Filter adata to only keep cells present in ref_barcodes
adata_filtered = adata[adata.obs_names.isin(ref_barcodes)].copy()

print(f"Original adata: {adata.n_obs} cells")
print(f"Filtered adata: {adata_filtered.n_obs} cells")

out_path = object_dir / f'RAW_slideseq_{adata_filtered.n_obs}_same_as_umap.h5ad'
adata_filtered.write(out_path)
