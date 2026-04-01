import scanpy as sc
import pandas as pd
import os

h5ad_path = "/scratch/mfafouti/Mommybrain/CAST/data/reference/Zeng/imputed/C57BL6J-638850-imputed-log2.h5ad"
out_dir = "/scratch/mfafouti/Mommybrain/CAST/data/reference/Zeng/imputed"

print(f"Loading {h5ad_path} ...")
adata = sc.read_h5ad(h5ad_path)
print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Save obs
obs_path = os.path.join(out_dir, "C57BL6J-638850-imputed-log2_obs.csv")
print(f"Saving obs to {obs_path} ...")
adata.obs.to_csv(obs_path)
print("obs saved.")

# Save each obsm slot
if len(adata.obsm) == 0:
    print("No obsm slots found.")
else:
    for key in adata.obsm.keys():
        obsm_path = os.path.join(out_dir, f"C57BL6J-638850-imputed-log2_obsm_{key}.csv")
        print(f"Saving obsm['{key}'] to {obsm_path} ...")
        pd.DataFrame(adata.obsm[key], index=adata.obs_names).to_csv(obsm_path)
        print(f"obsm['{key}'] saved.")

print("Done.")