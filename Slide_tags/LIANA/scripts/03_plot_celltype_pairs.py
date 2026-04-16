import argparse
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
import pandas as pd
import plotnine as p9
import liana as li
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot cell type pair interactions from lr_res')
parser.add_argument('--p_val', type=float, default=0.1,
                    help='interaction_padj threshold for significance filter')
parser.add_argument('--run_dir', type=str, required=True,
                    help='Path to the run output directory (e.g. out/runs/my_run)')
args = parser.parse_args()

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
out_dir = Path(args.run_dir) / 'liana_edgeR' / 'plots'
out_dir.mkdir(exist_ok=True, parents=True)

# =================== INPUT ===================
lr_res_path = Path(args.run_dir) / 'liana_edgeR' / 'lr_res.csv'
lr_res = pd.read_csv(lr_res_path)
print(lr_res.head())

# =================== FILTER ===================
significant_interactions = lr_res[lr_res['interaction_padj'] < args.p_val]
communication_counts = (
    significant_interactions
    .groupby(['source', 'target', 'interaction', 'interaction_stat'])
    .size()
    .reset_index(name='interaction_count')
    .sort_values(by='interaction_stat', ascending=False)
)
print(f"Significant interactions (padj < {args.p_val}): {len(significant_interactions)}")
print(communication_counts.head())

# =================== HEATMAP ===================
heatmap_data = communication_counts.pivot_table(
    index='source', columns='target', values='interaction_stat', aggfunc='mean'
)

plt.figure(figsize=(8, 6))
sns.heatmap(
    heatmap_data,
    cmap='RdBu_r',
    center=0,
    linewidths=0.5,
    linecolor='black',
    cbar_kws={'label': 'Log2FC (High = UP in CORT)'}
)
plt.title(f'CCC aggregated by cell type: interactions with padj<{args.p_val} shown')
plt.xlabel('Target Cell Type')
plt.ylabel('Source Cell Type')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()

heatmap_path = out_dir / f'interaction_heatmap_padj{args.p_val}.png'
plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
print(f"Heatmap saved to {heatmap_path}")
plt.close()

# =================== DOTPLOT ===================
max_abs_stat = significant_interactions['interaction_stat'].abs().max()

dotplot_fig = li.pl.dotplot(
    liana_res=significant_interactions,
    colour='interaction_stat',
    size='interaction_padj',
    inverse_size=True,
    orderby='interaction_stat',
    orderby_ascending=False,
    orderby_absolute=True,
    top_n=10,
    size_range=(0.5, 6),
    cmap='RdBu_r',
    figure_size=(19, 6)
)

dotplot_fig = dotplot_fig + p9.labs(x="", y="", title="") + \
    p9.scale_color_gradientn(
        colors=['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
                '#f7f7f7',
                '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f'],
        values=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        limits=(-max_abs_stat, max_abs_stat),
        name='Log2FC'
    )

dotplot_path = out_dir / f'dotplot_top10_padj{args.p_val}.png'
dotplot_fig.save(dotplot_path, dpi=300)
print(f"Dotplot saved to {dotplot_path}")
