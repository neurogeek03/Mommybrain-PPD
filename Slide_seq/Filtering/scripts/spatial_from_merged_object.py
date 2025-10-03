import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import distinctipy 
import plotly.express as px 
import plotly.graph_objects as go

# paths
adata_dir= '/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects'
output_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/spatial_plots_300_10_in_2'
os.makedirs(output_dir, exist_ok=True)

adata_path = os.path.join(adata_dir, "RAW_adata_filtered_220626_10_in_any_2_samples.h5ad")
type_column = "RCTD_first_type_rat"
spot_class_column = "RCTD_spot_class_rat"

# functions
def get_num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except:
        return -1

# Extract number before first "_" for ordering
def extract_num(celltype):
    try:
        return int(str(celltype).split("_", 1)[0])
    except ValueError:
        return -1  # fallback if no number

# build color palette
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)
color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)
color_df = color_df.sort_values("num_prefix")
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

# read data 
adata = sc.read_h5ad(adata_path)

# loop over all samples in the merged file 
for sample in adata.obs["sample"].unique():
    # subset to that sample
    adata_sample = adata[adata.obs["sample"] == sample]
    print(f"Processing sample: {sample}")
    # check if columns are there
    # check if columns are there
    singlets = adata_sample[adata_sample.obs[spot_class_column] == "singlet"].copy()
    if singlets.shape[0] == 0:
        print(f"⚠️ Skipping {sample}: no singlet spots")
        continue

    if type_column not in singlets.obs.columns or "X_spatial" not in singlets.obsm:
        print(f"⚠️ Skipping {sample}: missing data")
        continue
    # get coordinates 
    coords = singlets.obsm["X_spatial"]
    df = pd.DataFrame(coords, columns=["x", "y"], index=singlets.obs_names)
    second_type_df = singlets.obs[[type_column]].copy()
    merged_df = df.join(second_type_df, how='left')
    # Count cell types and identify top 30
    top30_celltypes = (
        merged_df[type_column]
        .value_counts()
        .nlargest(30)
        .index
    )
    # Work with string version to avoid category errors
    merged_df["celltype_plot"] = merged_df[type_column].astype(str)
    merged_df.loc[~merged_df["celltype_plot"].isin(top30_celltypes), "celltype_plot"] = "Other"
    # Set as categorical with desired order
    merged_df["celltype_plot"] = pd.Categorical(
        merged_df["celltype_plot"],
        categories=list(top30_celltypes) + ["Other"],
        ordered=True
    )
    # Sort categories based on number (ascending), with "Other" last
    ordered_celltypes = sorted(
        [ct for ct in merged_df["celltype_plot"].cat.categories if ct != "Other"],
        key=extract_num
    ) + ["Other"]
    # Apply new category order
    merged_df["celltype_plot"] = pd.Categorical(
        merged_df["celltype_plot"],
        categories=ordered_celltypes,
        ordered=True
    )
    # Build palette using CSV
    palette = {}
    for ct in ordered_celltypes:
        if ct == "Other":
            palette[ct] = "lightgrey"
        else:
            # fallback if celltype not found in CSV
            palette[ct] = label_to_hex.get(ct, "#000000")

    # ========== INDIVIDUAL PLOT ==========
    fig = px.scatter(
    merged_df,
    x="x",
    y="y",
    color="celltype_plot",
    color_discrete_map=palette,  # same palette dict you used
    hover_data={
            "RCTD_first_type_rat": True,
            "x": False,   # explicitly hide
            "y": False,   # explicitly hide
            "celltype_plot": False   # explicitly hide
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
        title=f"{sample} ({len(merged_df)} singlets, top 30 + Other)",
        xaxis_title="x",
        yaxis_title="y",
        legend_title=type_column, 
        # autosize = True, 
        width = 1500, 
        height = 1200
    )
    fig.update_layout(paper_bgcolor='white', plot_bgcolor='white', legend=go.layout.Legend(itemsizing='constant'))


    fig_path = os.path.join(output_dir, f"{type_column}_{sample}_filtered_spatial_RCTD.html")
    fig.write_html(fig_path)


    print(f"✅ Saved interactive plot to {fig_path}")

    # ============== PREVIOUS PNG CODE ==============

    # plt.figure(figsize=(12, 8))
    # ax = plt.gca()
    # sns.scatterplot(
    #     ax=ax,
    #     x="x",
    #     y="y",
    #     hue="celltype_plot",
    #     data=merged_df,
    #     s=4,
    #     alpha=0.8,
    #     linewidth=0,
    #     palette=palette
    # )


    # ax.set_title(f"{sample} ({len(merged_df)} singlets, top 30 + Other)", fontsize=14)
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
    # output_png = os.path.join(output_dir, f"{type_column}_{sample}_filtered_spatial_RCTD.png")
    # plt.savefig(output_png, dpi=300)
    # plt.close()
    # print(f"✅ Saved individual plot (top 30 + Other) to {output_png}")

