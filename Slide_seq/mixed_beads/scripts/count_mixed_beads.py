import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import scanpy as sc
import plotly.express as px

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

data_path = Path("/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/OUT/FINAL_RCTD_newgenelist/anndata_objects")

# =================== PARAMS ===================
sample = 'B01'

def classify_cell_type(cell_type):
    if isinstance(cell_type, str):
        ending = cell_type.split('_')[-1]
        if ending.endswith('NN'):
            return 'non-neuronal'
        elif ending in ['Glut', 'GABA']:
            return 'neuron'
    return 'other'

# =================== INPUT ===================
found_files = []
for file in data_path.iterdir():
    parts = file.name.split('_', 1)
    file_prefix = parts[0]
    if file_prefix == sample:
        found_files.append(file)
        print(f'Found matching file {file.stem}')

if found_files:
    adata = sc.read_h5ad(found_files[0])
    
    # Filter for beads with RCTD_spot_class => only the ones that were classified with RCTD
    adata = adata[adata.obs['RCTD_spot_class'].notna()].copy()

    # Classify cell types
    adata.obs['1st_type'] = adata.obs['RCTD_first_type_rat'].apply(classify_cell_type)
    adata.obs['2nd_type'] = adata.obs['RCTD_second_type_rat'].apply(classify_cell_type)
    
    # Count mixed rows
    mixed_beads_indices = (
        ((adata.obs['1st_type'] == 'neuron') & (adata.obs['2nd_type'] == 'non-neuronal')) |
        ((adata.obs['1st_type'] == 'non-neuronal') & (adata.obs['2nd_type'] == 'neuron'))
    )
    mixed_beads_count = mixed_beads_indices.sum()

    # Add a new metadata column for mixed beads
    adata.obs['mixed_bead_type'] = np.nan
    adata.obs.loc[mixed_beads_indices, 'mixed_bead_type'] = 'mixed'
    

# =================== Bead counts ===================
    print(f'Number of mixed beads: {mixed_beads_count}')
    all_rctd_beads_count = adata.obs['RCTD_spot_class'].notna().sum()
    print(f'Number of total RCTD beads: {all_rctd_beads_count}')
    ratio = (mixed_beads_count/all_rctd_beads_count)*100
    print(f'The mixed beads comprise {ratio:.2f}% of the sample')

# # =================== Spatial plot of mixed beads ===================
#     fig, ax = plt.subplots(figsize=(10, 10))
#     sc.pl.spatial(adata, color='mixed_bead_type', na_color='lightgrey',
#                   title=f"Spatial Plot for Sample {sample} - {ratio:.2f}% Mixed Beads",
#                   show=False, ax=ax, spot_size=100)
#     plot_file = output_base / f"{sample}_mixed_beads_spatial_plot.png"
#     plt.savefig(plot_file)
#     print(f"Saved spatial plot to {plot_file}")

# =================== Alluvial plot of first and second type ===================
    import plotly.graph_objects as go                                                                  
    import random   
    # Prepare data for alluvial plot
    alluvial_data = adata.obs.value_counts(['RCTD_first_type_rat', 'RCTD_second_type_rat']).reset_index(name='count')
    # Filter for the top 100 pairs
    alluvial_data = alluvial_data.sort_values(by='count', ascending=False).head(100)


    # # Debugging: Inspect the alluvial_data DataFrame
    # print("--- Alluvial Data Head ---")
    # print(alluvial_data.head())
    # print("\n--- Count Column Summary ---")
    # print(alluvial_data['count'].describe())

    # # Create the parallel categories plot
    # fig = px.parallel_categories(alluvial_data, dimensions=['RCTD_first_type_rat', 'RCTD_second_type_rat'],
    #                              color="RCTD_first_type_rat", color_continuous_scale=px.colors.sequential.Inferno,
    #                              title=f"Alluvial Plot for Sample {sample} (top 100 pairs)")

    # alluvial_plot_file = output_base / f"exp_{sample}_alluvial_plot.html"
    # fig.write_html(alluvial_plot_file)
    # print(f"Saved alluvial plot to {alluvial_plot_file}")

    # =================== GO SNAKEY ===================
    # Create nodes and links for Sankey diagram
    labels = list(pd.unique(alluvial_data[['RCTD_first_type_rat', 'RCTD_second_type_rat']].values.ravel('K')))
    source = [labels.index(s) for s in alluvial_data['RCTD_first_type_rat']]
    target = [labels.index(t) for t in alluvial_data['RCTD_second_type_rat']]
    value = alluvial_data['count']

    # Generate random colors for each label
    colors = ['#'+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(labels))]

    # Create the Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color=colors
        ),
        link=dict(
            source=source,
            target=target,
            value=value
        ))])

    fig.update_layout(title_text=f"Alluvial Plot for Sample {sample} (Top 100 Pairs)", font_size=10)
    alluvial_plot_file = output_base / f"100_{sample}_alluvial_plot.html"
    fig.write_html(alluvial_plot_file, include_plotlyjs='dist')
    print(f"Saved alluvial plot to {alluvial_plot_file}")