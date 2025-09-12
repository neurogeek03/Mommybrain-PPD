import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd

filename = 'umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad'
adata_path = f"/scratch/mfafouti/Mommybrain/Slide_seq/Integration/{filename}"
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Integration/umaps_filtered"        
os.makedirs(out_dir, exist_ok=True)

sample = 'all'

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


# Load harmonized AnnData
adata_all = sc.read_h5ad(adata_path)

# # FILTERING SUBCLASSES
# counts = adata_all.obs.groupby(['sample', 'RCTD_first_type_mouse']).size().unstack(fill_value=0)

# # Keep only cell types with at least 10 cells in **every** sample
# valid_celltypes = counts.columns[(counts >= 10).all(axis=0)]

# # Sum across samples to get overall prevalence and take top 40
# #top40 = counts[valid_celltypes].sum(axis=0).sort_values(ascending=False).head(40).index
# top40 = counts[valid_celltypes].sum(axis=0).sort_values(ascending=False).index

# # Subset the AnnData object
# adata_top40 = adata_all[adata_all.obs['RCTD_first_type_mouse'].isin(top40)].copy()

# print(f"Original cells: {adata_all.n_obs}, Filtered cells: {adata_top40.n_obs}")
# adata_all.write("/scratch/mfafouti/Mommybrain/Slide_seq/Integration/filtered_subclass_umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad")
# top40_counts = counts[top40]
# top40_counts.to_csv("top40_RCTD_first_type_mouse_counts.csv")

# adata_all = adata_top40

adata = sc.read_h5ad("/scratch/mfafouti/Mommybrain/Slide_seq/Integration/umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad")

# Count cells per sample and RCTD_first_type_mouse
counts = adata.obs.groupby(['sample', 'RCTD_first_type_mouse']).size().reset_index(name='count')

# Pivot table to make types as rows and samples as columns
pivot = counts.pivot(index='RCTD_first_type_mouse', columns='sample', values='count').fillna(0)

# Keep only the types where all samples have >=10 cells
types_to_keep = pivot[(pivot >= 10).all(axis=1)].index.tolist()

print(F"Types to keep: {types_to_keep}")
print("Number of types to keep: ", len(types_to_keep))

# Subset AnnData
adata_filtered = adata[adata.obs['RCTD_first_type_mouse'].isin(types_to_keep)].copy()

print(f"Original cells: {adata.n_obs}, Filtered cells: {adata_filtered.n_obs}")

adata_filtered.write("/scratch/mfafouti/Mommybrain/Slide_seq/Integration/filtered_10_subclass_sample_umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad")

adata_filtered = adata_all



sc.pl.umap(adata_all, color='RCTD_first_type_mouse', palette=label_to_hex, show=True)
plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
plt.savefig(os.path.join(out_dir, f"celltype_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
            dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()

sc.pl.umap(adata_all, color='sample', show=True)
plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
plt.savefig(os.path.join(out_dir, f"sample_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
            dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()

sc.pl.umap(adata_all, color='pregnancy',show=True)
plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
plt.savefig(os.path.join(out_dir, f"pregnancy_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
            dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()

sc.pl.umap(adata_all, color='treatment', show=True)
plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
plt.savefig(os.path.join(out_dir, f"treatment_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
            dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()