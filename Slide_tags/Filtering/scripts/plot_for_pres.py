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
import argparse

# python -i plot_for_pres.py tags /scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad /scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/report --mode top_n --top_n 30 --flip

# ============ ARGUMENT PARSING ============
parser = argparse.ArgumentParser(description='Generate UMAP and barplot figures from an AnnData object.')
parser.add_argument('method', choices=['seq','seq_nn','tags'], default='tags')
parser.add_argument('ad_path', type=Path, help='Path to the input .h5ad file.')
parser.add_argument('out_dir', type=Path,
                    help='Output directory. Defaults to the same directory as the input file.')
parser.add_argument('--mode', choices=['original', 'top_n'], default='original',
                    help='Barplot mode: original (all types) or top_n (filter to top N by count)')
parser.add_argument('--top_n', type=int, default=20,
                    help='Number of top cell types to show (only used in top_n mode)')
parser.add_argument('--flip', action='store_true',
                    help='Flip axes: horizontal bars instead of vertical')
args = parser.parse_args()

ad_path = args.ad_path

# ============ PARAMS ============
cell_type_column = 'subclass_name'
if args.method == 'seq':
    cell_type_column = 'RCTD_first_type_rat'
if args.method == 'seq_nn':
    cell_type_column = 'allen_class'

# ============ PATHS ============
out_dir = args.out_dir if args.out_dir is not None else ad_path.parent

out_dir.mkdir(parents=True, exist_ok=True)

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
    if args.method == 'seq':
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

if args.method == 'seq':
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
if args.method == 'seq':
    adata.obs[cell_type_column] = pd.Categorical(
        adata.obs[cell_type_column], 
        categories=ordered_labels, 
        ordered=True
    )

# 2. Create the color palette list in the correct order
if args.method == 'seq_nn':
    seq_nn_colors = {
        "30 Astro-Epen": "#F4A261",
        "31 OPC-Oligo":  "#4CC9F0",
        "32 OEC":        "#F28482",
        "33 Vascular":   "#A8D672",
        "34 Immune":     "#E78AC3"
    }
    scanpy_colors = [seq_nn_colors.get(label, '#cccccc') for label in adata.obs[cell_type_column].cat.categories]
else:
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
    save=scanpy_png_name, 
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

    if args.mode == 'top_n':
        top_indices = counts_df.nlargest(args.top_n, 'count').index
        counts_df = counts_df.loc[top_indices]

    if args.method == 'seq_nn':
        # --- SIMPLE MODE: sort by count, no grouping ---
        sorted_df = counts_df.sort_values('count', ascending=True).copy()
        x_positions = list(range(len(sorted_df)))
        x_labels = sorted_df['subclass_name'].tolist()
        sorted_df['x_pos'] = x_positions
    else:
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

        group_gap = 1.5

        for i, row in sorted_df.iterrows():
            current_group = row['group']

            if last_group is not None and current_group != last_group:
                current_pos += group_gap

            x_positions.append(current_pos)

            modified_name = str(row['subclass_name'])
            modified_name = re.sub(r'Glut', '', modified_name, flags=re.IGNORECASE)
            modified_name = re.sub(r'Gaba', '', modified_name, flags=re.IGNORECASE)
            modified_name = re.sub(r'NN', '', modified_name, flags=re.IGNORECASE)
            modified_name = re.sub(r'[\s_]+', ' ', modified_name)
            modified_name = modified_name.strip().strip('_')

            x_labels.append(modified_name)
            current_pos += 1
            last_group = current_group

        sorted_df['x_pos'] = x_positions

    # --- END: SORTING LOGIC ---

    if sorted_df.empty:
        print(f"Warning: No data to plot for '{barplot_column}'. Skipping barplot.")
    else:
        if args.method == 'seq_nn':
            seq_nn_colors = {
                "30 Astro-Epen": "#F4A261",
                "31 OPC-Oligo":  "#4CC9F0",
                "32 OEC":        "#F28482",
                "33 Vascular":   "#A8D672",
                "34 Immune":     "#E78AC3"
            }
            bar_colors = [seq_nn_colors.get(subcls, '#cccccc') for subcls in sorted_df['subclass_name']]
        else:
            bar_colors = [label_to_hex.get(subcls, '#cccccc') for subcls in sorted_df['subclass_name']]

        if args.flip:
            fig, ax = plt.subplots(figsize=(13, 2))
            ax.barh(
                sorted_df['x_pos'],
                sorted_df['count'],
                color=bar_colors
            )
            ax.set_xscale('log')
            ax.set_xlabel('Number of Cells (log scale)', fontsize=14)
            ax.set_yticks(sorted_df['x_pos'])
            ax.set_yticklabels(x_labels, fontsize=12)
            ax.set_ylim(min(x_positions) - 0.5, max(x_positions) + 0.5)
            ax.tick_params(axis='y', length=0)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
        else:
            fig, ax = plt.subplots(figsize=(20, 10))
            ax.bar(
                sorted_df['x_pos'],
                sorted_df['count'],
                color=bar_colors
            )
            ax.set_yscale('log')
            ax.set_ylabel('Number of Cells (log scale)', fontsize=14)
            ax.tick_params(axis='y', labelsize=14)
            ax.set_xticks(sorted_df['x_pos'])
            ax.set_xticklabels(x_labels, fontsize=12, rotation=90)
            ax.set_xlim(min(x_positions) - 0.5, max(x_positions) + 0.5)
            ax.tick_params(axis='x', length=0)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)

        ax.grid(False)
        plt.tight_layout()

        suffix = ''
        if args.mode == 'top_n':
            suffix += f'_top{args.top_n}'
        if args.flip:
            suffix += '_flipped'
        out_name = f'new_barplot_cells_per_subclass_log_vertical_grouped_spaced_min10_asc{suffix}.png'
        barplot_save_path = out_dir / out_name
        plt.savefig(barplot_save_path, dpi=300)
        plt.close()
        print(f"Saved barplot to {barplot_save_path}")

        # # ============ COMBINED FIGURE: UMAP + Barplot ============
        # if args.flip:
        #     combined_fig, axes = plt.subplots(
        #         2, 1, figsize=(15, 20),
        #         gridspec_kw={'height_ratios': [2, 1]}
        #     )
        # else:
        #     combined_fig, axes = plt.subplots(
        #         1, 2, figsize=(35, 12),
        #         gridspec_kw={'width_ratios': [2, 1]}
        #     )

        # # Draw UMAP directly into axes[0]
        # sc.pl.umap(
        #     adata,
        #     color=cell_type_column,
        #     title='',
        #     show=False,
        #     frameon=False,
        #     legend_loc=None,
        #     ax=axes[0]
        # )

        # # Redraw barplot into axes[1]
        # ax2 = axes[1]
        # if args.flip:
        #     ax2.barh(sorted_df['x_pos'], sorted_df['count'], color=bar_colors)
        #     ax2.set_xscale('log')
        #     ax2.set_xlabel('Number of Cells (log scale)', fontsize=14)
        #     ax2.set_yticks(sorted_df['x_pos'])
        #     ax2.set_yticklabels(x_labels, fontsize=12)
        #     ax2.set_ylim(min(x_positions) - 0.5, max(x_positions) + 0.5)
        #     ax2.tick_params(axis='y', length=0)
        #     ax2.spines['top'].set_visible(False)
        #     ax2.spines['right'].set_visible(False)
        #     ax2.spines['left'].set_visible(False)
        # else:
        #     ax2.bar(sorted_df['x_pos'], sorted_df['count'], color=bar_colors)
        #     ax2.set_yscale('log')
        #     ax2.set_ylabel('Number of Cells (log scale)', fontsize=14)
        #     ax2.tick_params(axis='y', labelsize=14)
        #     ax2.set_xticks(sorted_df['x_pos'])
        #     ax2.set_xticklabels(x_labels, fontsize=12, rotation=90)
        #     ax2.set_xlim(min(x_positions) - 0.5, max(x_positions) + 0.5)
        #     ax2.tick_params(axis='x', length=0)
        #     ax2.spines['top'].set_visible(False)
        #     ax2.spines['right'].set_visible(False)
        #     ax2.spines['bottom'].set_visible(False)
        # ax2.grid(False)

        # plt.tight_layout()
        # combined_save_path = out_dir / f'combined_umap_barplot{suffix}.png'
        # plt.savefig(combined_save_path, dpi=200, bbox_inches='tight')
        # plt.close()
        # print(f"Saved combined figure to {combined_save_path}")

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

# python -i plot_for_pres.py seq /scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/out/neurons/269646_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad /scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/report/seq --mode top_n --top_n 50 --flip