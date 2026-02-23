import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.patches import Patch
import plotly.express as px 
import plotly.graph_objects as go
import numpy as np
from pathlib import Path
import re

# ============ PARAMS ============
cell_type_column = 'RCTD_first_type_rat'

# ============ PATHS ============
project_folder = Path('/scratch/mfafouti/Mommybrain/Slide_seq/Filtering')
# in_dir = project_folder / 'NEW_list_merged_filtered'
# out_dir = project_folder / 'PRES'

out_dir = Path('/scratch/mfafouti/Mommybrain/Slide_seq/merging_NN_neurons/output/dim_red')
# in_dir = project_folder / 'out'

out_dir.mkdir(parents=True, exist_ok=True)

# ad_path = out_dir / "10_Celltype_subclass_223660_10_in_any_2_samples_UMAP.png"

# ad_path = Path('/scratch/mfafouti/Mommybrain/Slide_seq/Integration/FINAL_run_newgenelist/objects/adata_filtered_220626_10_in_any_2_samples_singlet_score_300.h5ad')

ad_path = Path('/scratch/mfafouti/Mommybrain/Slide_seq/merging_NN_neurons/output/dim_red/258708_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad')
# /scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/DE_after_mt_filter_108123_mincells_10_in_2_samples_slide_tags.h5ad

adata = sc.read_h5ad(ad_path)
# Replace 'your_column_name' with the actual name of the column you're interested in
column_name = cell_type_column

# Access the column in adata.obs and count unique values
unique_count = adata.obs[column_name].nunique()

# Print the result
print(f"Number of unique values in '{column_name}': {unique_count}")
file_stem = ad_path.stem
fig_path = out_dir / f"{cell_type_column}_{file_stem}.html"

# ============ GET COLORS ============
color_df_path = "/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv"
if not os.path.exists(color_df_path):
    print(f"Warning: Color file not found at {color_df_path}. Using default colors.")
    label_to_hex = {}
else:
    color_df = pd.read_csv(color_df_path, usecols=["name", "color_hex_triplet"])
    color_df["name"] = color_df["name"].str.replace(r"[ /-]", "_", regex=True)

    # Create mapping from label number (prefix before _) to hex color
    def get_num_prefix(label):
        try:
            return int(str(label).split("_")[0])
        except:
            return -1

    color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)

    # Sort CSV by numeric prefix
    color_df = color_df.sort_values("num_prefix")

    # Build dictionary mapping label to hex
    label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

# -------------- INTERACTIVE UMAP (PLOTLY) --------------
# (This section was already correct and is kept as-is)
umap_df = adata.obs.copy()
print(umap_df.head())
umap_df["UMAP1"] = adata.obsm["X_umap"][:, 0]
umap_df["UMAP2"] = adata.obsm["X_umap"][:, 1]
umap_df["label_number"] = umap_df[cell_type_column].str.split("_").str[0].astype(int)
ordered_labels = umap_df.groupby(cell_type_column, observed=True)["label_number"].first().sort_values().index.tolist()
umap_df[cell_type_column] = pd.Categorical(umap_df[cell_type_column], categories=ordered_labels, ordered=True)

umap_df[cell_type_column] = (umap_df[cell_type_column])

umap_df = umap_df.sort_values(cell_type_column)
umap_df.head()

# fig = px.scatter(
#         umap_df,
#         x="UMAP1",
#         y="UMAP2",
#         color=cell_type_column,
#         color_discrete_map=label_to_hex,  # <-- use your palette
#         hover_data={
#             cell_type_column: True,
#             "subclass_name": True,
#             "class_name": True,
#             "class_bootstrapping_probability": True,
#             "sample": True,
#             "treatment": True,
#             "UMAP1": False,   # explicitly hide
#             "UMAP2": False    # explicitly hide
#         }
#     )
# fig.update_traces(marker=dict(size=2, opacity=0.8))
# fig.update_layout(title=f"UMAP Plot Slide_tags data, n={adata.n_obs} colored by {cell_type_column}",
#                   paper_bgcolor='white', 
#                   plot_bgcolor='white', 
#                   legend=go.layout.Legend(itemsizing='constant'))
# fig.write_html(fig_path)

# print(f"Saved plotly umap to {fig_path}")


# ============ NEW SECTION: SCANPY UMAP PLOT ============
print("Generating Scanpy UMAP plot...")

# 1. Apply the same categorical ordering to the main adata.obs object
#    (Using the 'ordered_labels' list from the plotly section)
adata.obs[cell_type_column] = pd.Categorical(
    adata.obs[cell_type_column], 
    categories=ordered_labels, 
    ordered=True
)

# 2. Create the color palette list in the correct order
scanpy_colors = [label_to_hex.get(label, '#cccccc') for label in adata.obs[cell_type_column].cat.categories]

# 3. Assign the colors to adata.uns so scanpy can find them
adata.uns[f'{cell_type_column}_colors'] = scanpy_colors

# 4. Set the scanpy figure directory
sc.settings.figdir = str(out_dir)

# 5. Plot using sc.pl.umap
scanpy_fig_name = f"_{file_stem}_scanpy.svg"
sc.pl.umap(
    adata, 
    color=cell_type_column, 
    title=f"UMAP projection, n={adata.n_obs} cells", 
    save=scanpy_fig_name, 
    show=False,  # Don't show inline, just save
    frameon=False, # Common aesthetic choice
    legend_loc=None # <-- THIS PLACES LABELS ON CLUSTERS, REMOVING THE LEGEND
)

print(f"Saved scanpy umap to {out_dir / f'umap{scanpy_fig_name}'}")

scanpy_png_name = f"_{file_stem}_scanpy.png"
sc.pl.umap(
    adata, 
    color=cell_type_column, 
    title=f"UMAP projection, n={adata.n_obs} cells", 
    save=scanpy_fig_name, 
    show=False,  # Don't show inline, just save
    frameon=False, # Common aesthetic choice
    legend_loc=None # <-- THIS PLACES LABELS ON CLUSTERS, REMOVING THE LEGEND
)

print(f"Saved scanpy umap to {out_dir / f'umap{scanpy_png_name}'}")


# ============ Barplot of cells per subclass (VERTICAL LOG SCALE) ============
# ============ Barplot of cells per subclass (GROUPED, VERTICAL, LOG SCALE) ============
import re # <-- Added this import

barplot_column = cell_type_column
if barplot_column not in adata.obs.columns:
    print(f"Warning: '{barplot_column}' not in adata.obs. Skipping barplot.")
else:
    subclass_counts = adata.obs[barplot_column].value_counts()
    
    # --- START: SORTING LOGIC ---
    
    # 1. Convert Series to DataFrame for easier manipulation
    counts_df = subclass_counts.reset_index()
    counts_df.columns = ['subclass_name', 'count']
    counts_df = counts_df[counts_df['count'] > 10].copy()

    # 2. Define a function to extract the group
    def get_group(name):
        name_str = str(name).lower()
        if 'glut' in name_str:
            return 'Glut'
        elif 'gaba' in name_str:
            return 'Gaba'
        elif 'nn' in name_str or 'non-neuronal' in name_str:
            return 'NN'
        else:
            return 'Other'

    # 3. Create the 'group' column
    counts_df['group'] = counts_df['subclass_name'].apply(get_group)

    # 4. Filter to keep only the 3 specified groups
    filtered_df = counts_df[counts_df['group'].isin(['Glut', 'Gaba', 'NN'])].copy()

    # 5. Define the desired order of groups and apply it
    group_order = ['Glut', 'Gaba', 'NN']
    filtered_df['group'] = pd.Categorical(filtered_df['group'], categories=group_order, ordered=True)

    # 6. Sort by group (respecting categorical order), then by count (ascending)
    sorted_df = filtered_df.sort_values(by=['group', 'count'], ascending=[True, True])

    # Create new x-positions to introduce gaps
    x_positions = []
    x_labels = []
    current_pos = 0
    last_group = None
    
    # Define the size of the gap between groups
    group_gap = 1.5 # You can adjust this value (1.0 = normal bar width)
    
    for i, row in sorted_df.iterrows():
        current_group = row['group']
        
        # If the group changes, add the gap
        if last_group is not None and current_group != last_group:
            current_pos += group_gap
            
        x_positions.append(current_pos)
        
        # --- START: Label Modification (This is the change) ---
        # Remove the group names from the label string
        modified_name = str(row['subclass_name'])
        modified_name = re.sub(r'Glut', '', modified_name, flags=re.IGNORECASE)
        modified_name = re.sub(r'Gaba', '', modified_name, flags=re.IGNORECASE)
        modified_name = re.sub(r'NN', '', modified_name, flags=re.IGNORECASE)
        
        # Clean up leftover spaces or underscores
        modified_name = re.sub(r'[\s_]+', ' ', modified_name) # Replace 1+ spaces/underscores with one space
        modified_name = modified_name.strip() # Remove leading/trailing spaces
        modified_name = modified_name.strip('_') # Remove leading/trailing underscores
        
        x_labels.append(modified_name)
        # --- END: Label Modification ---
        
        current_pos += 1 # Advance position for the next bar
        last_group = current_group

    # Add the new positions to the DataFrame
    sorted_df['x_pos'] = x_positions

    # --- END: SORTING LOGIC ---

    if sorted_df.empty:
        print(f"Warning: No cells found for 'Glut', 'Gaba', or 'NN' in '{barplot_column}'. Skipping barplot.")
    else:
        # Adjusted figsize
        plt.figure(figsize=(20, 10)) 
        
        # --- Use the x_pos for plotting ---
        plt.bar(
            sorted_df['x_pos'], # <-- Corrected this
            sorted_df['count'],
            color=[label_to_hex.get(subcls, '#cccccc') for subcls in sorted_df['subclass_name']]
        )
        
        plt.yscale('log')
        plt.ylabel('Number of Cells (log scale)', fontsize=14) 
        
        plt.yticks(fontsize=14)
        
        # --- Use x_pos and cleaned x_labels for ticks ---
        plt.xticks(
            ticks=sorted_df['x_pos'], # <-- Corrected this
            labels=x_labels,          # <-- Corrected this
            fontsize=12, 
            rotation=90
        )
        
        plt.grid(False)

        # --- Correct xlim to use new positions ---
        plt.xlim(min(x_positions) - 0.5, max(x_positions) + 0.5)
        
        # --- Configure spines ---
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        # (Optional) Remove the x-axis ticks for a cleaner look
        ax.tick_params(axis='x', length=0)
        
        plt.tight_layout()
        
        # Updated file name
        barplot_save_path = out_dir / 'new_barplot_cells_per_subclass_log_vertical_grouped_spaced_min10_asc.png'
        plt.savefig(barplot_save_path, dpi=300)
        plt.close()
        print(f"Saved barplot to {barplot_save_path}")

# The print("Script finished.") line will be after this block

# The print("Script finished.") line will be after this block
    
# # ============ FIXED: Barplot of cells per subclass ============
# print("Generating bar plot...")

# # Use the column name from your original code. 
# # Note: This is 'MapMyCells_subclass_name', not 'subclass_name'
# barplot_column = cell_type_column
# if barplot_column not in adata.obs.columns:
#     print(f"Warning: '{barplot_column}' not in adata.obs. Skipping barplot.")
# else:
#     # <-- FIXED: Removed undefined 'my_order'. 
#     subclass_counts = adata.obs[barplot_column].value_counts()
    
#     # Sort by count (ascending) as in your original code
#     subclass_counts_sorted = subclass_counts.sort_values(ascending=True)
    
#     plt.figure(figsize=(10, 20))
    
#     # <-- FIXED: Replaced undefined 'color_dict' with 'label_to_hex'
#     plt.barh(
#         subclass_counts_sorted.index, 
#         subclass_counts_sorted.values,
#         color=[label_to_hex.get(subcls, '#cccccc') for subcls in subclass_counts_sorted.index]
#     )
#     plt.xscale('log')
#     plt.xlabel('Number of Cells')
    
#     # <-- FIXED: Removed undefined 'n_cells' from title
#     plt.title(f'Cell Number per Subclass')
    
#     plt.grid(False)  # turns off the grid

#     plt.xticks(fontsize=14)
#     plt.yticks(fontsize=10)

#     plt.gca().invert_yaxis() # Keep this if you want largest bar at the top

#     # 1. Remove top and bottom whitespace
#     # Set limits tight to the first and last bar
#     plt.ylim(len(subclass_counts_sorted) - 0.5, -0.5)

#     # 2. Remove plot borders (spines)
#     ax = plt.gca()
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.spines['left'].set_visible(False)
#     # The bottom spine ('x-axis') remains visible by default

#     plt.tight_layout()
    
#     barplot_save_path = out_dir / 'new_barplot_cells_per_subclass.svg'
#     plt.savefig(barplot_save_path, dpi=300)
#     plt.close()
#     print(f"Saved barplot to {barplot_save_path}")

# print("Script finished.")