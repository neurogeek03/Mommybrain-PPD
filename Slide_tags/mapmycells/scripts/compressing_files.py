import math
import anndata
import os
import scipy.sparse
import scanpy as sc 

input_dir = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender"

sample_list = ['BC3', 'BC9', 'BC15', 'BC14', 'BC28']

for sample in sample_list: 
    print(f"Processing {sample}...")

    sample_folder = os.path.join(input_dir, sample)
    h5ad_path = os.path.join(sample_folder, f'collapsed_mouse_id_{sample}.h5ad')

    ad = sc.read_h5ad(h5ad_path)
    count_matrix = ad.X
    obs = ad.obs
    var = ad.var

    file_size_bytes = os.path.getsize(h5ad_path)
    print("File size in bytes:", file_size_bytes)

    N = math.ceil(file_size_bytes/2147483648)
    num_rows = ad.shape[0]
    rows_per_subset = num_rows // N

    for i in range(N):
     start_idx = i * rows_per_subset
     end_idx = min((i + 1) * rows_per_subset, num_rows)

     ad = anndata.AnnData(
		X=count_matrix[start_idx:end_idx, :],
		obs=obs[start_idx:end_idx], #obs.iloc[start_idx:end_idx]
		var=var
	 )
     ad.write_h5ad(f'/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/mapmycells/{sample}_{i}_mapmycells.h5ad', compression='gzip')