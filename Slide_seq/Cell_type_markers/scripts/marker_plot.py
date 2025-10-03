import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.colors import LinearSegmentedColormap
import os 
import re

# ================== LOAD DATA ==================
path_base = "/scratch/mfafouti/Mommybrain/Slide_seq/Cell_type_markers/out/updated_marker_plot"
csv_path = os.path.join(path_base, "dotplot_data.csv")
spec_path = os.path.join(path_base, "subclass_gene_specificity.csv")
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Cell_type_markers/out/updated_marker_plot"
os.makedirs(out_dir, exist_ok=True) 

df = pd.read_csv(csv_path)
specificity = pd.read_csv(spec_path, index_col=0)

# Function to extract number from the string
def extract_number(name):
    match = re.search(r'\d+', name)
    return int(match.group()) if match else float('inf')  # put non-matching at end

# ================== GENE SELECTION ==================
# Define order of subclasses (if not predefined)
# nn_subclasses = df[df['group'].str.contains("Glut")]['group'].unique().tolist()  # or set manually if you want a specific order
# my_order = sorted(nn_subclasses, key=extract_number)
# print(my_order)

top_n = 2
ordered_genes = []
for subclass in df["group"].unique():
    if subclass in specificity.columns:
        top_genes = specificity[subclass].sort_values(ascending=False).head(top_n).index.tolist()
        ordered_genes.extend(top_genes)
    else:
        print(f"⚠️ Subclass '{subclass}' not in specificity data")

# Remove duplicates while preserving order
genes_for_plot = list(OrderedDict.fromkeys(ordered_genes))

# Restrict df to selected genes only
df = df[df["gene"].isin(genes_for_plot)].copy()

# Add general category of subclasses
df["category"] = df["group"].str.split("_").str[2]

# # Set categorical order for plotting
# df["gene"] = pd.Categorical(df["gene"], categories=genes_for_plot, ordered=True)
# df["group"] = pd.Categorical(df["group"], categories=my_order, ordered=True)
# # Sort for plotting
# df = df.sort_values(["gene", "group"])

# Sort genes for plotting
df["gene"] = pd.Categorical(df["gene"], categories=genes_for_plot, ordered=True)

print("Genes used for plot:", genes_for_plot)

# ================== CUSTOM PURPLE COLORMAP ==================
base_color = "#732B8B"
purple_palette = sns.light_palette(base_color, n_colors=6, input="hex", reverse=False)
purple_cmap = LinearSegmentedColormap.from_list("custom_purple", purple_palette)

# ===== Global Plot seetings =====
sns.set_theme(style="white")

# # ================== DOTPLOT ==================
# plt.figureure(figuresize=(30, 27))

# scatter = sns.scatterplot(
#     data=df,
#     x="gene",
#     y="group",
#     size="pct_expressing",
#     hue="mean_expression",
#     palette=purple_cmap,
#     sizes=(0, 200),
#     edgecolor="gray",
#     linewidth=0.5,
# )

# plt.xticks(rotation=90, fontsize=14)
# plt.yticks(fontsize=14)                   # subclass names
# plt.xlabel("Rattus Norvegicus Gene Symbol", fontsize=16)
# plt.ylabel("Allen Brain Institute Subclass", fontsize=16)
# plt.title("Dotplot of Marker Gene Expression", fontsize=18)
# plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
# plt.tight_layout(pad=4) 

# plt.savefigure(os.path.join(out_dir, f"top_{top_n}_genescustom_dotplot_from_csv.png"), dpi=300)

# ================== DOTPLOT (FACETED)==================
sns.set_theme(style="white")

g = sns.relplot(
    data=df,
    x="gene",
    y="group",
    size="pct_expressing",
    hue="mean_expression",
    palette=purple_cmap,
    sizes=(0, 200),
    edgecolor="gray",
    linewidth=0.5,
    col="category",
    col_order=["NN", "Glut", "Gaba"],
    facet_kws={"sharex": True, "sharey": False},  # 🔑 share x only
    height=9,
    aspect=1.2
)

g.set_xticklabels(rotation=90, fontsize=12)
g.set_ylabels("Allen Brain Institute Subclass", fontsize=14)
g.set_xlabels("Rattus Norvegicus Gene Symbol", fontsize=14)
g.set_titles("{col_name} cells")  # panel titles
g.figure.subplots_adjust(top=0.9)
g.figure.suptitle("Dotplot of Marker Gene Expression", fontsize=18)

plt.savefig(os.path.join(out_dir, f"top_{top_n}_genes_faceted_dotplot.png"), dpi=300)