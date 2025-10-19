import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
from pathlib import Path
import numpy as np

# ============ PARAMS ============
qc_column = 'class_bootstrapping_probability'

# ============ PATHS ============
project_folder = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Filtering')
in_dir = project_folder / 'NEW_list_merged_filtered'

ad_path = in_dir / "umap_filtered_150725_slide_tags.h5ad"

# adata = sc.read_h5ad(ad_path)

# Define the variables
x_col = "class_name"
y_col = "class_bootstrapping_probability"

file_stem = ad_path.stem
fig_path = in_dir / f"color_box_QC_{x_col}_{y_col}_{file_stem}.png"
csv_path = in_dir / f"QC_{x_col}_{y_col}_{file_stem}.csv"

hist_path = in_dir / f"Hist_{y_col}_{file_stem}.png"
# ============ Figure df ============
# df = adata.obs[[x_col, y_col]].copy()
# df.to_csv(csv_path)

df = pd.read_csv(csv_path)

# Order by median for readability (optional)
order = df.groupby(x_col)[y_col].median().sort_values(ascending=False).index

plt.figure(figsize=(10,6))
ax = sns.boxplot(
    x=x_col, 
    y=y_col, 
    data=df,
    order=order,
    palette="Set2",
    showcaps=True,
    boxprops={'facecolor':'none', 'linewidth':1.5},  # transparent box, colored edges
    medianprops={'color':'black', 'linewidth':1.5},
    whiskerprops={'linewidth':1.2},
    capprops={'linewidth':1.2},
    fliersize=0
)

# Make each box outline use the facecolor from the palette
for patch, color in zip(ax.patches, sns.color_palette("Set2", n_colors=len(order))):
    patch.set_edgecolor(color)

plt.xticks(rotation=45, ha='right')
plt.xlabel(x_col)
plt.ylabel(y_col)
plt.title(f"{y_col} by {x_col}")
plt.tight_layout()
plt.savefig(fig_path)



plt.figure(figsize=(6,4))
sns.histplot(x=y_col, data=df, bins=30, kde=False, color="steelblue")
plt.xlabel(y_col)
plt.ylabel("Count")
plt.title(f"Distribution of {y_col}")
plt.tight_layout()
plt.savefig(hist_path)
