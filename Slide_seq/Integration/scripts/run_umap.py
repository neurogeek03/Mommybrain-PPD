import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse
# ========== CONFIG ==========
# Paths
# filename = 'All_RCTD_types_singlet_score_0_slide_seq_15.h5ad'
# adata_path = f"/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects/{filename}"
# out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/OUT"        
# os.makedirs(out_dir, exist_ok=True)

# sample = filename.split('_')[1]

# ========== ARGUMENT PARSING ==========
parser = argparse.ArgumentParser(description='Process Slide-seq AnnData objects.')
parser.add_argument('input', type=str, help='Path to h5ad file')
parser.add_argument('outdir', type=str, help='Output directory')
args = parser.parse_args()

# Assign paths
adata_path = args.input
out_dir = args.outdir

# Create output directory and set Scanpy settings
os.makedirs(out_dir, exist_ok=True)
sc.settings.figdir = out_dir

meta_column = "RCTD_spot_class_rat"
filter_column = "RCTD_singlet_score_rat"
filter_value = 0

# sample = 'all_data'

# Load harmonized AnnData
adata_all = sc.read_h5ad(adata_path)
print(f'Number of cells: {adata_all.n_obs}')
n_reject = (adata_all.obs[meta_column] == "reject").sum()
print("Number of reject cells:", n_reject)
adata_all = adata_all[adata_all.obs[meta_column] != "reject"].copy()
n_unique = adata_all.obs[meta_column].nunique()
print("Number of unique values:", n_unique)
print(f'Number of cells: {adata_all.n_obs}')
#adata_all = adata_all[adata_all.obs[filter_column]>filter_value].copy()

# print(adata_all.obs["RCTD_second_type"].unique().tolist())
print('starting computations')
sc.pp.filter_cells(adata_all, min_counts=100)   # sensible thresholds
sc.pp.filter_genes(adata_all, min_cells=3)
sc.pp.normalize_total(adata_all, target_sum=1e4)
sc.pp.log1p(adata_all)
sc.pp.highly_variable_genes(adata_all, flavor='seurat_v3', n_top_genes=2000)
adata = adata_all[:, adata_all.var['highly_variable']]


print('done with log normalization')
sc.pp.scale(adata_all, max_value=10)
print('scaling done')
# sc.tl.pca(adata_all, svd_solver='arpack')
sc.tl.pca(adata_all, n_comps=50, svd_solver='randomized') 
print('pca done')
sc.pp.neighbors(adata_all, n_neighbors=10, n_pcs=30)
print('neighbors done')
#sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, method='pynndescent', n_jobs=16)
sc.tl.umap(adata_all)
print('done with umap!')
adata_all.write(os.path.join(out_dir, f"{adata_all.n_obs}_umap_filtered_{filter_value}_NEW_genelist_slide_seq_15.h5ad"))

# GET COLORS
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)

# Create mapping from label number (prefix before _) to hex color
def get_num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except:
        return -1

color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)

# Sort CSV by numeric prefix
color_df = color_df.sort_values("num_prefix")

# Build dictionary mapping label to hex
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

# ========== PLOTTING ==========

# 1. RCTD_second_type
sc.pl.umap(adata_all, color='RCTD_first_type_rat', palette=label_to_hex, show=True)
plt.title(f"UMAP of singlets after RCTD rat", fontsize=14)
plt.savefig(os.path.join(out_dir, f"all_filtered_{filter_value}_singlets_UMAP.png"),
            dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()

