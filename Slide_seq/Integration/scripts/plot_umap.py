import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.patches import Patch
import plotly.express as px 
import plotly.graph_objects as go
import numpy as np
import scipy.sparse as sp

filename = '1054147_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad' #CHANGE
# /scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects/
adata_path = f"/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects/{filename}"
object_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects"
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/OUT"        
os.makedirs(out_dir, exist_ok=True)

sample = 'filtered'
color_col = 'RCTD_singlet_score_rat'

# ============ GET COLORS ============
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

# ============ LOADING & FILTERING DATA ============
# Load harmonized AnnData
adata = sc.read_h5ad(adata_path)
n_unique = adata.obs["RCTD_spot_class_rat"].nunique()
print("Number of unique values:", n_unique)
# adata = adata[adata.obs["RCTD_singlet_score_rat"] > 300].copy()

# Count cells per sample and RCTD_first_type_rat
counts = adata.obs.groupby(['sample', 'RCTD_first_type_rat']).size().reset_index(name='count')
# Pivot table to make types as rows and samples as columns
pivot = counts.pivot(index='RCTD_first_type_rat', columns='sample', values='count').fillna(0)

# ----------------- Filter per cell type per sample -----------------
# # Keep only the types where all samples have >=10 cells
# types_to_keep = pivot[(pivot >= 10).all(axis=1)].index.tolist()

# # IGNOROING BAD SAMPLES 
# ignore_samples = ["B02", "B21"]
# samples_for_filtering = pivot.columns.difference(ignore_samples)
# types_to_keep = pivot[samples_for_filtering][(pivot[samples_for_filtering] >= 10).all(axis=1)].index.tolist()

types_to_keep = pivot[pivot.columns][(pivot[pivot.columns] >= 10).sum(axis=1) >= 2].index.tolist()

print(F"Types to keep: {types_to_keep}")
print("Number of types to keep: ", len(types_to_keep))

# Subset AnnData
adata_filtered = adata[adata.obs['RCTD_first_type_rat'].isin(types_to_keep)].copy()

print(f"Original cells: {adata.n_obs}, Filtered cells: {adata_filtered.n_obs}")

adata_filtered.obs["coronal_section"] = adata_filtered.obs["coronal_section"].replace({
    "early": "rostral",
    "late": "caudal"
})

ad_filtered_path = os.path.join(object_dir, f"RAW_adata_filtered_{adata_filtered.n_obs}_10_in_any_2_samples.h5ad")
adata_filtered.write(ad_filtered_path)

# adata_filtered = adata
num_cells = adata_filtered.n_obs

# ------------ Building legend ------------
type_column = "RCTD_first_type_rat"
type_counts = adata_filtered.obs[type_column].value_counts()

# Sanity check
n_unique = adata_filtered.obs[type_column].nunique()
n_unique = len(type_counts)
print(f"Unique {type_column}: {n_unique}")

top_types = type_counts.head(30).index.tolist()
# Build legend elements using your label_to_hex
legend_elements = [Patch(facecolor=label_to_hex.get(t, "#808080"), label=t) for t in top_types]
sc.pl.umap(adata_filtered, color='RCTD_first_type_rat', palette=label_to_hex, show=True)
plt.title(f"UMAP of singlets (n = {num_cells}) after RCTD mouse", fontsize=14)
# Build custom legend for top 30 only
plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5), title="Subclasses present with min n_cells=10 at all samples", fontsize=7)
plt.savefig(os.path.join(out_dir, f"FILTERED_celltype_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
            dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.close()

# -------------- INTERACTIVE UMAP --------------
umap_df = adata_filtered.obs.copy()
umap_df["UMAP1"] = adata_filtered.obsm["X_umap"][:, 0]
umap_df["UMAP2"] = adata_filtered.obsm["X_umap"][:, 1]
umap_df["label_number"] = umap_df["RCTD_first_type_rat"].str.split("_").str[0].astype(int)
ordered_labels = umap_df.groupby("RCTD_first_type_rat", observed=True)["label_number"].first().sort_values().index.tolist()
umap_df["RCTD_first_type_rat"] = pd.Categorical(umap_df["RCTD_first_type_rat"], categories=ordered_labels, ordered=True)

umap_df = umap_df.sort_values("RCTD_first_type_rat")
umap_df.head()

fig = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        color=color_col,
        color_discrete_map=label_to_hex,  # <-- use your palette
        hover_data={
            "RCTD_first_type_rat": True,
            "RCTD_second_type_rat": True,
            "RCTD_spot_class_rat": True,
            "RCTD_singlet_score_rat": True,
            "sample": True,
            "pregnancy": True,
            "treatment": True,
            "day": True,
            "coronal_section": True,
            "UMAP1": False,   # explicitly hide
            "UMAP2": False    # explicitly hide
        }
    )
fig.update_traces(marker=dict(size=2, opacity=0.8))
fig.update_layout(paper_bgcolor='white', plot_bgcolor='white', legend=go.layout.Legend(itemsizing='constant'))
fig_path = os.path.join(out_dir, f"{color_col}_BrainCanada_Slide-seq_umap_interactive.html")
fig.write_html(fig_path)


# sc.pl.umap(adata_filtered, color='RCTD_first_type_rat', palette=label_to_hex, show=True)
# plt.title(f"UMAP of singlets (n = {adata.n_obs}) after RCTD mouse", fontsize=14)
# plt.savefig(os.path.join(out_dir, f"20_Celltype_subclass_{adata.n_obs}_10_in_any_2_samples_UMAP.png"),
#             dpi=300, bbox_inches='tight', pad_inches=0.1)

# sc.pl.umap(adata, color='sample', show=True)
# plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
# plt.savefig(os.path.join(out_dir, f"sample_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
#             dpi=300, bbox_inches='tight', pad_inches=0.1)
# plt.close()

# sc.pl.umap(adata, color='pregnancy',show=True)
# plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
# plt.savefig(os.path.join(out_dir, f"pregnancy_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
#             dpi=300, bbox_inches='tight', pad_inches=0.1)
# plt.close()

# sc.pl.umap(adata, color='treatment', show=True)
# plt.title(f"UMAP of singlets after RCTD mouse", fontsize=14)
# plt.savefig(os.path.join(out_dir, f"treatment_subclass_{sample}_mouse_Slide-seq_umap_singlets.png"),
#             dpi=300, bbox_inches='tight', pad_inches=0.1)
# plt.close()