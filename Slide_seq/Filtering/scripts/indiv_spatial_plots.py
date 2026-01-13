import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
# import distinctipy 
import argparse
import plotly.express as px 
import plotly.graph_objects as go
from pandas.api.types import CategoricalDtype

# ========== ARGS ==========
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help="Path to input annotated anndata objects w/ RCTD output")
parser.add_argument("-o", "--output_dir", help="Path to output folder for the figures", default=None)
args = parser.parse_args()

filter_column = "RCTD_singlet_score_rat"
filter_value = 0
color = "RCTD_singlet_score_rat"

# ========== PATHS ==========
# adata_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/ratified_ref_anndata_objects"
adata_dir = args.input_dir
#output_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/spatial_plots_RCT_ratified_ref"
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

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

#print(label_to_hex)
print(adata_dir)
# ========== GET FILES ==========
h5ad_files = sorted([f for f in os.listdir(adata_dir) if f.endswith("_with_RCTD_mouse.h5ad")])
# h5ad_files = ["B37_with_RCTD_mouse.h5ad"]
spot_class_column = "RCTD_spot_class_rat"
type_column = "RCTD_first_type_rat"
print(f"Found {len(h5ad_files)} h5ad files")

for h5ad_file in h5ad_files:
    sample = h5ad_file.split("_")[0]
    print(f"📂 Reading: {sample}")

    adata = sc.read_h5ad(os.path.join(adata_dir, h5ad_file))

    #filtering
    adata = adata[adata.obs[filter_column]>filter_value].copy()

    # if spot_class_column not in adata.obs.columns:
    #     print(f"⚠️ Skipping {sample}: no RCTD_spot_class_rat")
    #     continue

    # singlets = adata[adata.obs[spot_class_column].notna()].copy()
    # Define a list of the categories you want to keep
    categories_to_keep = ["singlet", "doublet_certain", "doublet_uncertain"]

    # Use the .isin() method to filter for all categories in the list
    singlets = adata[adata.obs[spot_class_column].isin(categories_to_keep)].copy()

    # if singlets.shape[0] == 0:
    #     print(f"⚠️ Skipping {sample}: no singlet spots")
    #     continue

    # if type_column not in singlets.obs.columns or "X_spatial" not in singlets.obsm:
    #     print(f"⚠️ Skipping {sample}: missing data")
    #     continue

    coords = singlets.obsm["X_spatial"]
    df = pd.DataFrame(coords, columns=["x", "y"], index=singlets.obs_names)
    second_type_df = singlets.obs[[type_column, spot_class_column, color]].copy()
    
    merged_df = df.join(second_type_df, how='left')

    # Count cell types and identify top 30
    # top30_celltypes = (
    #     merged_df[type_column]
    #     .value_counts()
    #     .nlargest(30)
    #     .index
    # )

    # # Work with string version to avoid category errors
    # merged_df["celltype_plot"] = merged_df[type_column].astype(str)
    # merged_df.loc[~merged_df["celltype_plot"].isin(top30_celltypes), "celltype_plot"] = "Other"

    # # Set as categorical with desired order
    # merged_df["celltype_plot"] = pd.Categorical(
    #     merged_df["celltype_plot"],
    #     # categories=list(top30_celltypes) + ["Other"],
    #     ordered=True
    # )

    # Extract number before first "_" for ordering
    def extract_num(celltype):
        try:
            return int(str(celltype).split("_", 1)[0])
        except (ValueError, IndexError):
            return 99999
    
    unique_celltypes = merged_df[type_column].unique()

    # Sort the unique list to create your desired order
    ordered_celltypes = sorted(
        unique_celltypes,
        key=extract_num
    )
    # Create a new CategoricalDtype with your specific order
    cat_type = CategoricalDtype(categories=ordered_celltypes, ordered=True)

    # Update the column in the DataFrame to use this new ordered type
    # This tells pandas the "correct" way to sort this column
    merged_df[type_column] = merged_df[type_column].astype(cat_type)

    # --- 3. Sort the Entire DataFrame ---

    # Now, .sort_values() will use your custom order automatically
    merged_df_sorted = merged_df.sort_values(by=type_column)
    # # Sort categories based on number (ascending), with "Other" last
    # ordered_celltypes = sorted(
    #     [ct for ct in merged_df["celltype_plot"].cat.categories if ct != "Other"],
    #     key=extract_num
    # ) + ["Other"]

    # -----------------------
    # Build palette using CSV
    # -----------------------
    # palette = {}
    # for ct in ordered_celltypes:
    #     if ct == "Other":
    #         palette[ct] = "lightgrey"
    #     else:
    #         # fallback if celltype not found in CSV
    #         palette[ct] = label_to_hex.get(ct, "#000000")

    # ========== INDIVIDUAL PLOT ==========
    fig = px.scatter(
    merged_df_sorted,
    x="x",
    y="y",
    color=color,
    color_discrete_map=label_to_hex,  # same palette dict you used
    hover_data={
            "RCTD_first_type_rat": True,
            "RCTD_spot_class_rat": True,
            "x": False,   # explicitly hide
            "y": False,   # explicitly hide

        }
    )
    # match style
    fig.update_traces(
        marker=dict(size=4, opacity=0.8, line=dict(width=0))
    )

    # invert y-axis
    fig.update_yaxes(autorange="reversed")

    # titles and labels
    fig.update_layout(
        title=f"{sample} ({len(merged_df_sorted)} beads, colored by {color})",
        xaxis_title="x",
        yaxis_title="y",
        legend_title=type_column, 
        # autosize = True, 
        width = 1500, 
        height = 1200
    )
    fig.update_layout(paper_bgcolor='white', plot_bgcolor='white', legend=go.layout.Legend(itemsizing='constant'))


    fig_path = os.path.join(output_dir, f"full_no_rej_{color}_{sample}_filtered_spatial_RCTD.html")
    fig.write_html(fig_path)


    print(f"✅ Saved interactive plot to {fig_path}")

    # # ========== INDIVIDUAL PLOT ==========
    # plt.figure(figsize=(12, 8))
    # ax = plt.gca()

    # sns.scatterplot(
    #     ax=ax,
    #     x="x",
    #     y="y",
    #     hue=type_column,
    #     data=merged_df,
    #     s=4,
    #     alpha=0.8,
    #     linewidth=0,
    #     palette=label_to_hex
    # )

    # ax.set_title(f"{sample} ({len(merged_df)} singlets, all cell types)", fontsize=14)
    # ax.set_xlabel("x")
    # ax.set_ylabel("y")
    # ax.invert_yaxis()

    # # Custom legend with larger dots
    # handles, labels = ax.get_legend_handles_labels()
    # new_handles = [
    #     Line2D([0], [0], marker='o', color='w', label=label,
    #            markerfacecolor=handle.get_color(), markersize=8)
    #     for handle, label in zip(handles, labels)
    # ]
    # ax.legend(
    #     handles=new_handles,
    #     bbox_to_anchor=(1.05, 1),
    #     loc='upper left',
    #     borderaxespad=0,
    #     title=type_column
    # )

    # plt.tight_layout()
    # output_png = os.path.join(output_dir, f"{type_column}_ratified_coronal_ref_{sample}_spatial_RCTD.png")
    # plt.savefig(output_png, dpi=300)
    # plt.close()
    # print(f"✅ Saved individual plot (top 30 + Other) to {output_png}")


