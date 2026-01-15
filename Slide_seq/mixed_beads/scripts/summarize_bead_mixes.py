import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# =================== INPUT ===================
csv_path = project_path / 'out' / 'combined_run' /'382810_20260114_mixed_bead_counts_summary.csv'

# =================== PROCESSING ===================
df = pd.read_csv(csv_path)

# Prepare data for stacking
labels = df['sample']
neuron_to_non_neuronal = df['neuron_to_non-neuronal_count']
non_neuronal_to_neuron = df['non-neuronal_to_neuron_count']
total_beads = df['total_rctd_beads']
other_beads = total_beads - (neuron_to_non_neuronal + non_neuronal_to_neuron)

# Calculate percentages
neuron_to_non_neuronal_perc = (neuron_to_non_neuronal / total_beads) * 100
non_neuronal_to_neuron_perc = (non_neuronal_to_neuron / total_beads) * 100
other_beads_perc = (other_beads / total_beads) * 100

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8))

# Stacked bar
# Stacked bar
ax.bar(labels, neuron_to_non_neuronal_perc, label='Neuron to Non-Neuronal')
ax.bar(labels, non_neuronal_to_neuron_perc, bottom=neuron_to_non_neuronal_perc, label='Non-Neuronal to Neuron')
ax.bar(labels, other_beads_perc, bottom=neuron_to_non_neuronal_perc + non_neuronal_to_neuron_perc, label='Other', color='gray')


# Labels and title
ax.set_ylabel('Percentage (%)')
ax.set_xlabel('Sample')
ax.set_title('Mixed Bead Counts per Sample (Normalized)')
ax.legend()
plt.xticks(rotation=45)
plt.tight_layout()

# =================== OUTPUT ===================
output_plot_path = output_base / 'combined_run' / f'{stamp}_mixed_bead_counts_summary.png'
plt.savefig(output_plot_path)
print(f"Plot saved to {output_plot_path}")

print(f"------ Script finished at {datetime.now().strftime('%M%S%H_%Y%m%d')} ------")