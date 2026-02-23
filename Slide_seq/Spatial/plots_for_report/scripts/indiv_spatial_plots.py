import os
import re
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# ========== ARGS ==========
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_adata", help="Path to combined anndata object (all samples)")
parser.add_argument("-o", "--output_dir", help="Path to output folder for the figures", default=None)
parser.add_argument("-s", "--sample", help="Sample name to process (if omitted, all samples are processed)", default=None)
args = parser.parse_args()

# ========== PATHS ==========
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

# ========== GET COLORS ==========
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])
color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)

def get_num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except:
        return -1

color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)
color_df = color_df.sort_values("num_prefix")
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

# ========== BROAD CLASS CLASSIFICATION ==========
def normalize(name):
    return re.sub(r"[ /-]", "_", str(name))

def infer_broad_class(name):
    n = normalize(name)
    if n.endswith("_NN"):
        return "Non-neuronal"
    elif "_Glut" in n:
        return "Glutamatergic"
    elif "_IMN" in n:
        return "Glutamatergic"
    elif "_Gaba" in n:
        return "GABAergic"
    else:
        return "Unknown"

BROAD_CLASSES = ["Glutamatergic", "GABAergic", "Non-neuronal"]

# ========== LOAD COMBINED ADATA ==========
print(f"Reading combined adata: {args.input_adata}")
adata_all = sc.read_h5ad(args.input_adata)

spot_class_column = "RCTD_spot_class_rat"
type_column = "RCTD_first_type_rat"

# ========== DETERMINE SAMPLES TO PROCESS ==========
if args.sample is not None:
    samples = [args.sample]
else:
    samples = sorted(adata_all.obs["sample"].unique().tolist())

print(f"Processing {len(samples)} sample(s): {samples}")

# ========== HELPERS ==========
def extract_num(celltype):
    try:
        return int(str(celltype).split("_", 1)[0])
    except ValueError:
        return -1


def make_plot(ax, merged_df, ordered_celltypes, palette, title, type_column):
    sns.scatterplot(
        ax=ax,
        x="x",
        y="y",
        hue="celltype_plot",
        data=merged_df,
        s=6,
        alpha=0.8,
        linewidth=0,
        palette=palette
    )
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.get_legend().remove()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.margins(0)


# ========== PROCESS EACH SAMPLE ==========
for sample in samples:
    print(f"Processing: {sample}")

    adata = adata_all[adata_all.obs["sample"] == sample].copy()

    if spot_class_column not in adata.obs.columns:
        print(f"  Skipping {sample}: no {spot_class_column}")
        continue

    singlets = adata[adata.obs[spot_class_column] == "singlet"].copy()
    if singlets.shape[0] == 0:
        print(f"  Skipping {sample}: no singlet spots")
        continue

    if type_column not in singlets.obs.columns or "X_spatial" not in singlets.obsm:
        print(f"  Skipping {sample}: missing {type_column} or X_spatial")
        continue

    coords = singlets.obsm["X_spatial"]
    df = pd.DataFrame(coords, columns=["x", "y"], index=singlets.obs_names)
    merged_df = df.join(singlets.obs[[type_column]].copy(), how="left")
    merged_df["celltype_str"] = merged_df[type_column].astype(str)
    merged_df["broad_class"] = merged_df["celltype_str"].apply(infer_broad_class)

    # ---- Figure 1: all cells, top 30 + Other ----
    top30_celltypes = (
        merged_df[type_column]
        .value_counts()
        .nlargest(30)
        .index
    )

    merged_df["celltype_plot"] = merged_df["celltype_str"].copy()
    merged_df.loc[~merged_df["celltype_plot"].isin(top30_celltypes), "celltype_plot"] = "Other"

    ordered_celltypes = sorted(
        [ct for ct in top30_celltypes],
        key=extract_num
    ) + ["Other"]

    merged_df["celltype_plot"] = pd.Categorical(
        merged_df["celltype_plot"],
        categories=ordered_celltypes,
        ordered=True
    )

    palette = {ct: (label_to_hex.get(ct, "#000000") if ct != "Other" else "lightgrey")
               for ct in ordered_celltypes}

    fig, ax = plt.subplots(figsize=(5, 5))
    make_plot(ax, merged_df, ordered_celltypes, palette,
              f"{sample} ({len(merged_df)} singlets, top 30 + Other)", type_column)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.93)
    out_png = os.path.join(output_dir, f"{type_column}_ratified_coronal_ref_{sample}_spatial_RCTD.png")
    plt.savefig(out_png, dpi=300, bbox_inches=None)
    plt.close()
    print(f"  Saved: {out_png}")

    # ---- Figures 2-4: per broad class, colored cells on gray background ----
    for broad_class in BROAD_CLASSES:
        class_mask = merged_df["broad_class"] == broad_class
        n_colored = class_mask.sum()

        # All spots shown; only the target class gets color, rest are gray
        plot_df = merged_df.copy()
        plot_df["celltype_plot"] = plot_df["celltype_str"].copy()

        # Identify cell types in this broad class
        class_celltypes = sorted(
            plot_df.loc[class_mask, "celltype_plot"].unique().tolist(),
            key=extract_num
        )

        # Everything outside the class becomes "Background"
        plot_df.loc[~class_mask, "celltype_plot"] = "Background"

        ordered = class_celltypes + ["Background"]
        plot_df["celltype_plot"] = pd.Categorical(
            plot_df["celltype_plot"],
            categories=ordered,
            ordered=True
        )

        class_palette = {ct: label_to_hex.get(ct, "#000000") for ct in class_celltypes}
        class_palette["Background"] = "lightgrey"

        # Plot background first (gray), then colored on top
        fig, ax = plt.subplots(figsize=(10, 10))

        bg_df = plot_df[plot_df["celltype_plot"] == "Background"]
        ax.scatter(bg_df["x"], bg_df["y"], s=25, alpha=0.4, linewidths=0,
                   color="lightgrey", zorder=1)

        fg_df = plot_df[plot_df["celltype_plot"] != "Background"]
        for ct in class_celltypes:
            ct_df = fg_df[fg_df["celltype_plot"] == ct]
            if ct_df.empty:
                continue
            ax.scatter(ct_df["x"], ct_df["y"], s=25, alpha=0.8, linewidths=0,
                       color=class_palette[ct], zorder=2)

        ax.set_title(f"{sample} — {broad_class} ({n_colored} spots)", fontsize=14)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.margins(0)

        fig.subplots_adjust(left=0, right=1, bottom=0, top=0.93)
        class_tag = broad_class.replace("-", "_").replace(" ", "_")
        out_png = os.path.join(output_dir, f"{type_column}_{class_tag}_{sample}_spatial_RCTD.png")
        plt.savefig(out_png, dpi=300, bbox_inches=None)
        plt.close()
        print(f"  Saved: {out_png}")
