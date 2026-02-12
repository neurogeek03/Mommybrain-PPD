import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Define paths
project_path = Path.cwd().parents[0]
data_path = project_path / 'out' / 'celltype_interact' / 'comm_counts_005_interaction_stat.csv'
output_dir = project_path / 'out' / 'figures'
output_dir.mkdir(exist_ok=True, parents=True)
output_heatmap_path = output_dir / 'interaction_heatmap.png'

# Load the data
try:
    df = pd.read_csv(data_path)
except FileNotFoundError:
    print(f"Error: The file {data_path} was not found.")
    exit()

# Pivot the table to create a matrix for the heatmap
# Rows: source, Columns: target, Values: interaction_stat
# If there are multiple interactions between the same source-target pair,
# we'll take the mean of their interaction_stat for simplicity in this visualization.
# You might want to adjust this aggregation method if a different one is preferred (e.g., max, sum).
heatmap_data = df.pivot_table(index='source', columns='target', values='interaction_stat', aggfunc='mean')

# Create the heatmap
plt.figure(figsize=(12, 10)) # Adjust figure size as needed
sns.heatmap(
    heatmap_data,
    cmap='RdBu_r', # Red-Blue reversed colormap: Red for high positive, Blue for high negative
    center=0,      # Center the colormap at 0, useful for diverging stats
    linewidths=.5, # Add lines between cells for better separation
    linecolor='black',
    cbar_kws={'label': 'Interaction Statistic'} # Label for the color bar
)

plt.title('Interaction Statistic Heatmap by Cell Type Pair')
plt.xlabel('Target Cell Type')
plt.ylabel('Source Cell Type')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout() # Adjust layout to prevent labels from overlapping

# Save the heatmap
plt.savefig(output_heatmap_path, dpi=300, bbox_inches='tight')
print(f"Heatmap saved to {output_heatmap_path}")

# Display plot (optional, useful if running interactively)
# plt.show()
