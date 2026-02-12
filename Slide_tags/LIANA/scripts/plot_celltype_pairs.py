import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import seaborn as sns
import glob

import plotnine as p9
import liana as li

project_path = Path.cwd().parents[0]
lr_res_path = project_path / 'out'/ 'first_try' / 'lr_liana_res.csv'
output_dir = project_path / 'out' / 'celltype_interact'
output_dir.mkdir(exist_ok=True, parents=True)

lr_res = pd.read_csv(lr_res_path)

# pd.set_option('display.max_columns', None, index_col=0)
print(lr_res.head())

# Option 1: Count the number of significant interactions
p_val = 0.05
significant_interactions = lr_res[lr_res['interaction_padj'] < p_val]
communication_counts = significant_interactions.groupby(['source', 'target', 'interaction','interaction_stat']).size().reset_index(name='interaction_count')
communication_counts = communication_counts.sort_values(by='interaction_stat', ascending=False)
print("--- Communication Counts (significant interactions) ---")
print(communication_counts.head())

# Create the heatmap
output_heatmap_path = output_dir / f'interaction_heatmap_{p_val}_.png'

heatmap_data = communication_counts.pivot_table(index='source', columns='target', values='interaction_stat', aggfunc='mean')

# Create the heatmap
plt.figure(figsize=(8, 6)) # Adjust figure size as needed
sns.heatmap(
    heatmap_data,
    cmap='RdBu_r', # Red-Blue reversed colormap: Red for high positive, Blue for high negative
    center=0,      # Center the colormap at 0, useful for diverging stats
    linewidths=.5, # Add lines between cells for better separation
    linecolor='black',
    cbar_kws={'label': 'Interaction Statistic (High = UP in CORT)'} # Label for the color bar
)

plt.title(f'CCC aggregated by cell type: interactions with padj<{p_val} shown')
plt.xlabel('Target Cell Type')
plt.ylabel('Source Cell Type')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout() # Adjust layout to prevent labels from overlapping

# Save the heatmap
plt.savefig(output_heatmap_path, dpi=300, bbox_inches='tight')
print(f"Heatmap saved to {output_heatmap_path}")

# Visualizing cell types along with pathways 
print("\n--- Inspecting 'interaction_padj' ---")
print(significant_interactions['interaction_padj'].head())
print(f"Min 'interaction_padj': {significant_interactions['interaction_padj'].min()}")
print(f"Max 'interaction_padj': {significant_interactions['interaction_padj'].max()}")
significant_interactions.head()

dotplot_fig = li.pl.dotplot(liana_res=significant_interactions,
                     colour='interaction_stat',
                     size='interaction_padj',
                     inverse_size=True,
                     orderby='interaction_stat',
                     orderby_ascending=False,
                     orderby_absolute=True,
                     top_n=10,
                     size_range=(0.5, 6),
                     cmap='RdBu_r',
                     figure_size=(12, 5)
                     )

dotplot_fig = dotplot_fig + p9.labs(x="", y="", title="")

# # Option 2: Average the interaction_stat for significant interactions
# communication_scores = significant_interactions.groupby(['source', 'target'])['interaction_stat'].mean().reset_index(name='mean_interaction_stat')
# communication_scores = communication_scores.sort_values(by='mean_interaction_stat', ascending=False)
# print("\n--- Mean Interaction Score (significant interactions) ---")
# print(communication_scores.head())

dotplot_fig_path =  output_dir / f'p_{p_val}_dotsize_TOP10_dotplot_figure.png'
dotplot_fig.save(dotplot_fig_path, dpi=300)