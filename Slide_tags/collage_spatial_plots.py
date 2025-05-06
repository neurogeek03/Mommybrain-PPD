import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import anndata as ad
import pandas as pd

coords_dir = '/project/s/shreejoy/hudsonhu/RatAllCTOuts/spatial_coords_bender'
out_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Spatial'
project_path = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags"
adata_dir = os.path.join(project_path, "Post_bender")
csv_path = os.path.join(adata_dir, "hex_codes_subclass.csv")

sample_list = ["BC13", "BC14", "BC15", "BC28", "BC3", "BC9"]

# Load cell type color dict
df = pd.read_csv(csv_path, header=None, names=['cell_type', 'color'])
color_dict = dict(zip(df['cell_type'], df['color']))

# Create 2x3 subplot grid
fig, axes = plt.subplots(2, 3, figsize=(30, 12))

for i, sample in enumerate(sample_list):
    row = i // 3
    col = i % 3
    ax = axes[row, col]

    adata_file = os.path.join(adata_dir, f"{sample}", f"{sample}_spatial.h5ad")
    ad = sc.read_h5ad(adata_file)

    # Reading and subsetting coordinates
    coords_file = os.path.join(coords_dir, f'coords_{sample}.txt')
    coords = pd.read_csv(coords_file, header=0, delim_whitespace=True, quotechar='"')
    coords_subset = coords.iloc[:, :3].set_index(coords.columns[0])

    # Extract spatial data from obs
    obs_df = ad.obs.copy()
    obs_df['cell_bc'] = ad.obs_names
    valid_obs = obs_df[(obs_df['x_um'] != 0) | (obs_df['y_um'] != 0)]

    valid_coords = valid_obs[['x_um', 'y_um']]
    valid_cell_class = valid_obs['MapMyCells_cell_type']

    # Plot into subplot
    sns.scatterplot(
        ax=ax,
        x=valid_coords['x_um'],
        y=valid_coords['y_um'],
        hue=valid_cell_class,
        palette=color_dict,
        s=5,
        alpha=0.7,
        legend=False  # Suppress legends in subplots
    )

    ax.set_title(f'{sample} ({len(valid_coords)} cells)', fontsize=14)
    ax.set_xlabel('X Coordinate (µm)')
    ax.set_ylabel('Y Coordinate (µm)')

# Add a single legend for the whole figure (based on the last sample's classes)
handles, labels = ax.get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    title='Cell Class',
    bbox_to_anchor=(1.02, 0.5),
    loc='center left',
    borderaxespad=0.,
    markerscale=5
)

# Final layout and save
plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave room for the legend
plt.savefig(os.path.join(out_dir, "combined_spatial_all_samples.png"), dpi=300)
plt.show()
