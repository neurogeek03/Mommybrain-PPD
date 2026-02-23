import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.patches import Patch
import plotly.express as px 
import plotly.graph_objects as go
import numpy as np
from pathlib import Path

# ============ PARAMS ============
cell_type_column = 'RCTD_first_type_rat'
method = 'Slide-seq'

# ============ PATHS ============
project_folder = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Filtering')
# in_dir = project_folder / 'NEW_list_merged_filtered'
out_dir = project_folder / 'out'

# ad_path = out_dir / "PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
ad_path = '/scratch/mfafouti/Mommybrain/Slide_seq/merging_NN_neurons/output/dim_red/258708_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad'

adata = sc.read_h5ad(ad_path)
file_stem = ad_path.stem
fig_path = out_dir / f"{cell_type_column}_{file_stem}.html"

# ============ GET COLORS ============
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])
# color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)

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

# -------------- INTERACTIVE UMAP --------------
umap_df = adata.obs.copy()
print(umap_df.head())
umap_df["UMAP1"] = adata.obsm["X_umap"][:, 0]
umap_df["UMAP2"] = adata.obsm["X_umap"][:, 1]
umap_df["label_number"] = umap_df[cell_type_column].str.split(" ").str[0].astype(int)
ordered_labels = umap_df.groupby(cell_type_column, observed=True)["label_number"].first().sort_values().index.tolist()
umap_df[cell_type_column] = pd.Categorical(umap_df[cell_type_column], categories=ordered_labels, ordered=True)

umap_df[cell_type_column] = (umap_df[cell_type_column])

umap_df = umap_df.sort_values(cell_type_column)
umap_df.head()

fig = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        color=cell_type_column,
        color_discrete_map=label_to_hex,  # <-- use your palette
        hover_data={
            cell_type_column: True,
            "subclass_name": True,
            "class_name": True,
            "class_bootstrapping_probability": True,
            "sample": True,
            "treatment": True,
            "UMAP1": False,   # explicitly hide
            "UMAP2": False    # explicitly hide
        }
    )
fig.update_traces(marker=dict(size=2, opacity=0.8))
fig.update_layout(title=f"UMAP Plot {method} data, n={adata.n_obs} colored by {cell_type_column}",
                  paper_bgcolor='white', 
                  plot_bgcolor='white', 
                  legend=go.layout.Legend(itemsizing='constant'))
fig.write_html(fig_path)

print(f"Saved umap to {ad_path}")