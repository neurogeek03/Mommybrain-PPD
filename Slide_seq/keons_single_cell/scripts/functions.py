import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as pd
import scanpy as sc
import seaborn as sns

def classify_cell_type(cell_type):
    if isinstance(cell_type, str):
        ending = cell_type.split('_')[-1]
        if ending.endswith('NN'):
            return 'non-neuronal'
        elif ending in ['Glut', 'GABA']:
            return 'neuron'
    return 'other'

def plot_cell_counts(obj, obs_column, title, filename, use_log=True):
    sns.set_style("whitegrid")
    plt.figure(figsize=(50, 6))
    
    # 1. Get counts (already sorted high-to-low by default)
    counts = obj.obs[obs_column].value_counts()
    
    # 2. Create barplot
    ax = sns.barplot(
        x=counts.index, 
        y=counts.values, 
        order=counts.index, 
        palette='magma'
    )
    
    # 3. Conditional Log Scale
    if use_log:
        ax.set_yscale("log")
        ylabel = 'Number of Beads (Log10)'
        title_suffix = "(Log Scale)"
    else:
        ylabel = 'Number of Beads (Linear)'
        title_suffix = "(Linear Scale)"
    
    # 4. Formatting
    plt.title(f"{title} {title_suffix}")
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(ylabel)
    plt.xlabel('Cell Type')
    
    # Save to your defined directory
    plt.tight_layout()
    plt.savefig(f"{sc.settings.figdir}/{filename}", dpi=300)
