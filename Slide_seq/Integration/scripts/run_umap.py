import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
# ========== CONFIG ==========
# Paths
# filename ='class_B03_with_RCTD_mouse.h5ad'
# adata_path = f"/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/anndata_objects/{filename}" 
filename = 'B32_with_RCTD_mouse.h5ad'
adata_path = f"/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/anndata_objects/{filename}"
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Integration/umaps"        
os.makedirs(out_dir, exist_ok=True)

sample = filename.split('_')[1]

meta_column = "RCTD_spot_class_mouse"

# sample = 'all_data'

# Load harmonized AnnData
adata_all = sc.read_h5ad(adata_path)
adata_all = adata_all[adata_all.obs[meta_column] == "singlet"].copy()

# print(adata_all.obs["RCTD_second_type"].unique().tolist())
print('starting computations')
sc.pp.highly_variable_genes(adata_all, flavor='seurat_v3', n_top_genes=2000)
sc.pp.normalize_total(adata_all, target_sum=1e4)
sc.pp.log1p(adata_all)

print('done with log normalization')
sc.pp.scale(adata_all, max_value=10)
sc.tl.pca(adata_all, svd_solver='arpack')
sc.pp.neighbors(adata_all, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_all)
print('done with umap!')
# adata_all.write("/scratch/mfafouti/Mommybrain/Slide_seq/Integration/umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad")

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

# # Custom color dictionary for RCTD_second_type
# color_dict = {
#     "001_CLA_EPd_CTX_Car3_Glut": "#ffff00",
#     "004_L6_IT_CTX_Glut": "#00ffff",
#     "006_L4_5_IT_CTX_Glut": "#ff00ff",
#     "007_L2_3_IT_CTX_Glut": "#ff0000",
#     "009_L2_3_IT_PIR_ENTl_Glut": "#006666",
#     "012_MEA_Slc17a7_Glut": "#0080ff",
#     "014_LA_BLA_BMA_PA_Glut": "#660000",
#     "016_CA1_ProS_Glut": "#33cc33",
#     "017_CA3_Glut": "#663300",
#     "022_L5_ET_CTX_Glut": "#000099",
#     "025_CA2_FC_IG_Glut": "#00ffcc",
#     "030_L6_CT_CTX_Glut": "#cc9966",
#     "037_DG_Glut": "#003366",
#     "045_OB_STR_CTX_Inh_IMN": "#9999ff",
#     "061_STR_D1_Gaba": "#666666",
#     "062_STR_D2_Gaba": "#330000",
#     "093_RT_ZI_Gnb3_Gaba": "#666666",
#     "101_ZI_Pax6_Gaba": "#330000",
#     "151_TH_Prkcd_Grin2c_Glut": "#003300",
#     "152_RE_Xi_Nox4_Glut": "#ffd9e6",
#     "318_Astro_NT_NN": "#0000ff",
#     "319_Astro_TE_NN": "#666699",
#     "325_CHOR_NN": "#ff007f",
#     "326_OPC_NN": "#666666",
#     "327_Oligo_NN": "#990000",
#     "330_VLMC_NN": "#999900",
#     "331_Peri_NN": "#33cccc",
#     "333_Endo_NN": "#ffcc99",
#     "334_Microglia_NN": "#ffb6c1",
#     "338_Lymphoid_NN": "#cc33cc",
# }

# ========== PLOTTING ==========

# 1. RCTD_second_type
sc.pl.umap(adata_all, color='RCTD_first_type_mouse', palette=label_to_hex, show=True)
plt.title(f"UMAP of singlets after RCTD rat", fontsize=14)
plt.savefig(os.path.join(out_dir, f"subclas_{sample}_MOUSE_Slide-seq_umap_{sample}_singlets.png"),
            dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()

# # 2. Sample (autocolor)
# sc.pl.umap(adata_all, color='sample', show=False)
# plt.title("Harmony - color: sample ID", fontsize=14)
# plt.savefig(os.path.join(out_dir, "Slide-seq_umap_Harmony_sampleID.png"),
#             dpi=300, bbox_inches='tight', pad_inches=0.1)
# plt.close()

# # 3. Treatment (autocolor)
# sc.pl.umap(adata_all, color='treatment', show=False)
# plt.title("Harmony - color: treatment", fontsize=14)
# plt.savefig(os.path.join(out_dir, "Slide-seq_umap_Harmony_treatment.png"),
#             dpi=300, bbox_inches='tight', pad_inches=0.1)
# plt.close()

# # 4. Pregnant/nulliparous (autocolor)
# sc.pl.umap(adata_all, color='pregnant', show=False)
# plt.title("Harmony - color: pregnant/nulliparous", fontsize=14)
# plt.savefig(os.path.join(out_dir, "Slide-seq_umap_Harmony_pregnant.png"),
#             dpi=300, bbox_inches='tight', pad_inches=0.1)
# plt.close()