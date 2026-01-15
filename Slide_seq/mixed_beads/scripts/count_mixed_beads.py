import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import scanpy as sc
import plotly.express as px
import seaborn as sns

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out' / 'combined_run'
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
results = []
found_files = []
for file in data_path.iterdir():
    parts = file.name.split('_', 1)
    file_prefix = parts[0]
    # print(f'Found sample: {sample}')
    # found_files.append(file)
    if file_prefix == sample:
        found_files.append(file)
    print(f'Found matching file {file.stem}')

for sample_file in found_files:
    parts = sample_file.name.split('_', 1)
    sample = parts[0]
    adata = sc.read_h5ad(sample_file)
    
    # Filter for beads with RCTD_spot_class => only the ones that were classified with RCTD
    adata = adata[adata.obs['RCTD_spot_class'].notna()].copy()

    # Classify cell types
    adata.obs['1st_type'] = adata.obs['RCTD_first_type_rat'].apply(classify_cell_type)
    adata.obs['2nd_type'] = adata.obs['RCTD_second_type_rat'].apply(classify_cell_type)
    
    # Define and count different mixed bead types
    neuron_to_non_neuronal_indices = (adata.obs['1st_type'] == 'neuron') & (adata.obs['2nd_type'] == 'non-neuronal')
    non_neuronal_to_neuron_indices = (adata.obs['1st_type'] == 'non-neuronal') & (adata.obs['2nd_type'] == 'neuron')

    neuron_to_non_neuronal_count = neuron_to_non_neuronal_indices.sum()
    non_neuronal_to_neuron_count = non_neuronal_to_neuron_indices.sum()

    # Add a new metadata column for mixed beads
    adata.obs['mixed_bead_type'] = np.nan
    adata.obs.loc[neuron_to_non_neuronal_indices, 'mixed_bead_type'] = 'neuron_to_non-neuronal'
    adata.obs.loc[non_neuronal_to_neuron_indices, 'mixed_bead_type'] = 'non_neuronal_to_neuron'
    

# =================== Bead counts ===================
    # print(f'Number of mixed beads: {mixed_beads_count}')
    all_rctd_beads_count = adata.obs['RCTD_spot_class'].notna().sum()
    print(f'Number of mixed beads (neuron to non-neuronal): {neuron_to_non_neuronal_count}')
    ratio_n_nn = (neuron_to_non_neuronal_count / all_rctd_beads_count) * 100 if all_rctd_beads_count > 0 else 0
    print(f'The neuron to non-neuronal mixed beads comprise {ratio_n_nn:.2f}% of the sample')

    print(f'Number of mixed beads (non-neuronal to neuron): {non_neuronal_to_neuron_count}')
    ratio_nn_n = (non_neuronal_to_neuron_count / all_rctd_beads_count) * 100 if all_rctd_beads_count > 0 else 0
    print(f'The non-neuronal to neuron mixed beads comprise {ratio_nn_n:.2f}% of the sample')

    results.append({
        'sample': sample,
        'neuron_to_non-neuronal_count': neuron_to_non_neuronal_count,
        'non-neuronal_to_neuron_count': non_neuronal_to_neuron_count,
        'ratio_neuron_to_non-neuronal': ratio_n_nn,
        'ratio_non-neuronal_to_neuron': ratio_nn_n,
        'total_rctd_beads': all_rctd_beads_count
    })
# =================== Spatial plot of mixed beads ===================
    fig, ax = plt.subplots(figsize=(10, 10))
    sc.pl.spatial(adata, color='mixed_bead_type', na_color='lightgrey',
                  title=f"Spatial Plot for Sample {sample}",
                  show=False, ax=ax, spot_size=100)
    plot_file = output_base / f"{sample}_mixed_beads_spatial_plot.png"
    plt.savefig(plot_file)
    print(f"Saved spatial plot to {plot_file}")

    # =================== Heatmap of mixed beads ===================
    # Get the mixed beads
    mixed_beads_obs = adata.obs[(adata.obs['1st_type'] == 'neuron') & (adata.obs['2nd_type'] == 'non-neuronal') |
                                (adata.obs['1st_type'] == 'non-neuronal') & (adata.obs['2nd_type'] == 'neuron')]

    # Create a dataframe with neuron and non-neuronal types
    heatmap_data = []
    for _, row in mixed_beads_obs.iterrows():
        if row['1st_type'] == 'neuron':
            neuron_type = row['RCTD_first_type_rat']
            non_neuronal_type = row['RCTD_second_type_rat']
        else:
            neuron_type = row['RCTD_second_type_rat']
            non_neuronal_type = row['RCTD_first_type_rat']
        heatmap_data.append({'neuron': neuron_type, 'non_neuronal': non_neuronal_type})

    if heatmap_data:
        heatmap_df = pd.DataFrame(heatmap_data)
        
        # Create a contingency table
        contingency_table = pd.crosstab(heatmap_df['non_neuronal'], heatmap_df['neuron'])

        # Create the heatmap
        plt.figure(figsize=(12, 10))
        sns.heatmap(contingency_table, annot=True, fmt='d', cmap='viridis')
        plt.title(f'Heatmap of Neuron to Non-Neuronal Bead Mixes for Sample {sample}')
        plt.xlabel('Neuron Types')
        plt.ylabel('Non-Neuronal Types')
        
        heatmap_plot_file = output_base / f"{sample}_mixed_beads_heatmap.png"
        plt.savefig(heatmap_plot_file, bbox_inches='tight')
        plt.close()
        print(f"Saved heatmap to {heatmap_plot_file}")
# =================== Alluvial plot of first and second type ===================
    import plotly.graph_objects as go                                                                  
    import random   
    # Prepare data for alluvial plot
    alluvial_data = adata.obs.value_counts(['RCTD_first_type_rat', 'RCTD_second_type_rat']).reset_index(name='count')
    # Filter for the top 100 pairs
    alluvial_data = alluvial_data.sort_values(by='count', ascending=False).head(10)


    # # Debugging: Inspect the alluvial_data DataFrame
    # print("--- Alluvial Data Head ---")
    # print(alluvial_data.head())
    # print("\n--- Count Column Summary ---")
    # print(alluvial_data['count'].describe())

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
    alluvial_plot_file = output_base / f"10_{sample}_alluvial_plot.html"
    fig.write_html(alluvial_plot_file, include_plotlyjs='dist')
    print(f"Saved alluvial plot to {alluvial_plot_file}")

    # =================== HEATMAP ===================


# # saving summary csv for all samples 
# results_df = pd.DataFrame(results)
# csv_output_path = output_base / f"{stamp}_mixed_bead_counts_summary.csv"
# results_df.to_csv(csv_output_path, index=False)
# print(f"Saved summary of mixed bead counts to {csv_output_path}")

