import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from pathlib import Path


BROAD_CLASS_ORDER = ['Glutamatergic', 'GABAergic', 'Non-neuronal']


PAN_FIRST = ['Slc17a7', 'Gad2']

def build_gene_order(gene_panel_map):
    """
    Return gene list with pan markers (Slc17a7, Gad2) on rows 1-2, then the
    remaining genes ordered by panel (Glutamatergic → GABAergic → Non-neuronal).
    """
    gene_order = [g for g in PAN_FIRST if g in gene_panel_map]
    for bc in BROAD_CLASS_ORDER:
        for g, panel in gene_panel_map.items():
            if panel == bc and g not in gene_order:
                gene_order.append(g)
    return gene_order


def load_subclass_colors(color_csv=None, subclasses=None):
    """
    Load hex colors from a cluster_annotation CSV (columns: name, color_hex_triplet).
    Falls back to tab20 colormap if the file is absent or not provided.
    Returns {normalized_name: color}.
    """
    if color_csv is not None:
        p = Path(color_csv)
        if p.exists():
            df = pd.read_csv(p, usecols=['name', 'color_hex_triplet'])
            return dict(zip(df['name'], df['color_hex_triplet']))
        print(f"WARNING: color_csv not found: {p}. Using tab20 fallback.")

    cmap = plt.get_cmap('tab20')
    return {sub: cmap(i % 20) for i, sub in enumerate(subclasses or [])}


def draw_dotplot_panel(ax, gene_order, subclasses, mean_piv, pct_capped,
                       vmin, vmax, subclass_colours,
                       expr_cmap='Purples', dot_scale=200, bar_height=0.8):
    """
    Draw one dotplot panel onto ax.
    Returns the scatter object (needed for the shared colorbar).
    """
    n_genes = len(gene_order)
    bar_y = n_genes + 0.5

    xs, ys, sizes, colours = [], [], [], []
    for gi, gene in enumerate(gene_order):
        for si, sub in enumerate(subclasses):
            xs.append(si)
            ys.append(gi)
            sizes.append(float(pct_capped.at[gene, sub]) * dot_scale)
            colours.append(float(mean_piv.at[gene, sub]))

    sc = ax.scatter(xs, ys, s=sizes, c=colours, cmap=expr_cmap,
                    vmin=vmin, vmax=vmax, edgecolors='none')

    for si, sub in enumerate(subclasses):
        colour = subclass_colours.get(sub, '#808080')
        ax.add_patch(Rectangle(
            (si - 0.4, bar_y), 0.8, bar_height,
            facecolor=colour, edgecolor='none', clip_on=False,
        ))

    ax.set_xlim(-0.5, len(subclasses) - 0.5)
    ax.set_ylim(bar_y + bar_height + 0.2, -0.5)

    for gi in range(n_genes):
        ax.axhline(gi, color='#f0f0f0', linewidth=0.3, zorder=0)
    for si in range(len(subclasses) - 1):
        ax.axvline(si + 0.5, color='#d0d0d0', linewidth=0.4, zorder=0)

    return sc


def add_size_legend(ax, dot_scale, dot_max, legend_pcts=None, bbox=(1.12, 0.85)):
    if legend_pcts is None:
        legend_pcts = [0.2, 0.4, 0.6, 0.8]
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
               markersize=np.sqrt(p * dot_scale) * 0.6, label=f'{int(p * 100)}%')
        for p in legend_pcts if p <= dot_max
    ]
    ax.legend(handles=handles, title='% expressing', loc='upper left',
              bbox_to_anchor=bbox, frameon=False, fontsize=8, title_fontsize=9)
