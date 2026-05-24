import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
import pandas as pd
import seaborn as sns

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
out_dir = project_path / 'out' / 'liana_edgeR' / 'plots'
out_dir.mkdir(exist_ok=True, parents=True)

# =================== INPUT ===================
# Pre-aggregated communication counts from 03_plot_celltype_pairs.py or manual export
data_path = project_path / 'out' / 'celltype_interact' / 'comm_counts_005_interaction_stat.csv'
df = pd.read_csv(data_path)

# =================== HEATMAP ===================
# Rows: source, Columns: target, Values: mean interaction_stat across pairs
heatmap_data = df.pivot_table(
    index='source', columns='target', values='interaction_stat', aggfunc='mean'
)

plt.figure(figsize=(32, 10))
sns.heatmap(
    heatmap_data,
    cmap='RdBu_r',
    center=0,
    linewidths=0.5,
    linecolor='black',
    cbar_kws={'label': 'Interaction Statistic'}
)
plt.title('Interaction Statistic Heatmap by Cell Type Pair')
plt.xlabel('Target Cell Type')
plt.ylabel('Source Cell Type')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()

heatmap_path = out_dir / 'interaction_heatmap.png'
plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
print(f"Heatmap saved to {heatmap_path}")
plt.close()
