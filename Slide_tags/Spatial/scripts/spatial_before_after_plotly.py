#!/usr/bin/env python3
# ========== IMPORTS ==========
import scanpy as sc
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

# ========== PARAMETERS ==========
sample_list = ["BC13", "BC14", "BC28", "BC15", "BC3", "BC9"]
celltype_col = 'subclass_name'

# ========== PATHS ==========
mommybrain_folder = Path('/project/rrg-shreejoy/MommyBrain/Slide_tags/Pipeline_data/')
project_folder = Path('/scratch/mfafouti/Mommybrain/Slide_tags')
coords_dir = mommybrain_folder / 'spatial_coordinates'
conf_positioned_dir = project_folder / 'Spatial' / 'conf_positioned'
out_dir = project_folder / 'Spatial' / 'figures' / 'updated'
adata_dir = project_folder / 'Filtering' / 'out'
adata_path = adata_dir / "PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
out_dir.mkdir(exist_ok=True, parents=True)
html_path = out_dir / f"spatial_before_after_{celltype_col}.html"
fig_path_before = out_dir / f"before_filtering_{celltype_col}_combined_spatial.png"
fig_path_after = out_dir / f"conf_positioned_filtered_{celltype_col}_combined_spatial.png"

# ========== READING INTEGRATED ADATA FILE ==========
adata_all = sc.read_h5ad(adata_path)
print(adata_all.obs.head())
print(adata_all)

# ============ GET COLORS ============
color_df = pd.read_csv("/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv", usecols=["name", "color_hex_triplet"])

def get_num_prefix(label):
    try:
        return int(str(label).split("_")[0])
    except:
        return -1

color_df["num_prefix"] = color_df["name"].apply(get_num_prefix)
color_df = color_df.sort_values("num_prefix")
label_to_hex = dict(zip(color_df["name"], color_df["color_hex_triplet"]))

# ============ LOAD CONF_POSITIONED BARCODES PER SAMPLE ============
def load_conf_barcodes(sample):
    path = conf_positioned_dir / f"trekker_{sample}_{sample}_barcodes_ConfPositionedNuclei.tsv"
    if not path.exists():
        print(f"Warning: No conf_positioned file for {sample}")
        return None
    barcodes = pd.read_csv(path, header=None, names=['barcode'])
    barcodes['barcode'] = barcodes['barcode'].str.replace('-1', '', regex=False)
    return set(barcodes['barcode'])

# ============ PRE-COMPUTE DATA PER SAMPLE ============
sample_data = {}
for sample in sample_list:
    sample_ad = adata_all[adata_all.obs['sample'] == sample].copy()
    df = sample_ad.obs[[celltype_col]].reset_index(names='barcode').copy()
    df['barcode'] = df['barcode'].str.replace('-1', '', regex=False)

    coords_path = coords_dir / f'coords_{sample}.csv'
    if not coords_path.exists():
        print(f"Warning: Coords file not found for {sample}. Skipping.")
        sample_data[sample] = None
        continue

    obs_df = pd.read_csv(coords_path)
    merged_df = df.merge(obs_df, left_on='barcode', right_on='cell_bc', how='left')
    merged_df.set_index('barcode', inplace=True)

    # Before: all cells with valid spatial coordinates (same filter as existing script)
    before_df = merged_df[(merged_df['x_um'] != 0) | (merged_df['y_um'] != 0)].copy()
    before_df = before_df[~before_df[celltype_col].str.contains('NN', na=False)]

    # After: intersect with conf_positioned barcodes
    conf_barcodes = load_conf_barcodes(sample)
    if conf_barcodes is not None:
        after_df = before_df[before_df.index.isin(conf_barcodes)].copy()
    else:
        print(f"Warning: No conf_positioned barcodes for {sample}. After = Before.")
        after_df = before_df.copy()

    sample_data[sample] = {'before': before_df, 'after': after_df}
    print(f"{sample}: before={len(before_df)} cells, after={len(after_df)} cells")

# ============ BUILD SUBPLOT TITLES WITH CELL COUNTS ============
subplot_titles = []
for sample in sample_list:
    if sample_data[sample] is not None:
        n_before = len(sample_data[sample]['before'])
        n_after = len(sample_data[sample]['after'])
        subplot_titles.extend([
            f"{sample} — Before (n={n_before})",
            f"{sample} — After ConfPositioned (n={n_after})"
        ])
    else:
        subplot_titles.extend([
            f"{sample} — Before (no data)",
            f"{sample} — After (no data)"
        ])

# ============ CREATE FIGURE ============
fig = make_subplots(
    rows=len(sample_list),
    cols=2,
    subplot_titles=subplot_titles,
    horizontal_spacing=0.05,
    vertical_spacing=0.04
)

# Track legend entries — legendgroup ensures clicking toggles all traces of that cell type
legend_added = set()

for i, sample in enumerate(sample_list):
    row = i + 1
    if sample_data[sample] is None:
        continue

    for col, condition in enumerate(['before', 'after'], start=1):
        plot_df = sample_data[sample][condition]

        for cell_type, group in plot_df.groupby(celltype_col):
            color = label_to_hex.get(cell_type, '#cccccc')
            show_legend = cell_type not in legend_added
            if show_legend:
                legend_added.add(cell_type)

            fig.add_trace(
                go.Scatter(
                    x=group['x_um'],
                    y=group['y_um'],
                    mode='markers',
                    marker=dict(size=3, color=color, opacity=0.7),
                    name=cell_type,
                    legendgroup=cell_type,
                    showlegend=show_legend,
                    text=group[celltype_col],
                    hovertemplate='%{text}<br>x: %{x:.1f} µm<br>y: %{y:.1f} µm<extra></extra>'
                ),
                row=row, col=col
            )

# ============ LAYOUT ============
fig.update_layout(
    height=500 * len(sample_list),
    width=1200,
    title_text=f"Spatial plots — Before vs After ConfPositioned filtering ({celltype_col})",
    title_font_size=16,
    legend=dict(
        title=celltype_col,
        itemsizing='constant',
        font=dict(size=9),
        tracegroupgap=2
    ),
    plot_bgcolor='white',
    paper_bgcolor='white'
)

fig.update_xaxes(
    showgrid=False,
    zeroline=False,
    showline=True,
    linecolor='black',
    scaleanchor=None
)
fig.update_yaxes(
    showgrid=False,
    zeroline=False,
    showline=True,
    linecolor='black',
    scaleanchor=None
)

# Make each subplot square by matching y axis scale to x axis
for i, sample in enumerate(sample_list):
    for col in [1, 2]:
        fig.update_yaxes(scaleanchor=f'x{(i * 2 + col) if (i * 2 + col) > 1 else ""}',
                         scaleratio=1,
                         row=i + 1, col=col)

fig.write_html(html_path)
print(f"Saved interactive HTML to {html_path}")

# ============ STATIC MATPLOTLIB PLOTS ============
def save_static_plot(condition, fig_path):
    fig_static, axes = plt.subplots(2, 3, figsize=(20, 10))

    for i, sample in enumerate(sample_list):
        row = i // 3
        col = i % 3
        ax = axes[row, col]

        if sample_data[sample] is None:
            ax.set_title(f'{sample} (no data)', fontsize=14)
            ax.axis('off')
            continue

        plot_df = sample_data[sample][condition]

        sns.scatterplot(
            ax=ax,
            x=plot_df['x_um'],
            y=plot_df['y_um'],
            hue=plot_df[celltype_col],
            palette=label_to_hex,
            s=5,
            alpha=0.7,
            legend=False
        )

        ax.set_title(f'{sample} ({len(plot_df)} cells)', fontsize=14)
        ax.set_xlabel('X Coordinate (µm)')
        ax.set_ylabel('Y Coordinate (µm)')

    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    print(f"Saved static figure to {fig_path}")
    plt.close(fig_static)

save_static_plot('before', fig_path_before)
save_static_plot('after', fig_path_after)