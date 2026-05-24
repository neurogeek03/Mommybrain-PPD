"""
Step 2: Read marker_dotplot_data.csv and produce a 3-panel dotplot.

Panels: Glutamatergic | GABAergic | Non-neuronal
Gene order follows MARKER_CATEGORIES definition in utils/markers.py.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / 'code'))

from utils.markers import get_gene_panel_map
from utils.plot_utils import (
    BROAD_CLASS_ORDER,
    build_gene_order,
    load_subclass_colors,
    draw_dotplot_panel,
    add_size_legend,
)

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ 02_plot_markers started at {stamp} ------")

# =================== CONFIG ===================
CSV_PATH  = PROJECT / 'out' / 'marker_dotplot_data.csv'
COLOR_CSV = Path('/scratch/mfafouti/CAN-2026-poster/data/cluster_annotation_term.csv')
FIG_WIDTH   = 12
FIG_HEIGHT  = 14
EXPR_CMAP   = 'Purples'
DOT_MAX     = 0.8
DOT_SCALE   = 200
STD_SCALE   = True   # per-gene 0–1 scaling of mean expression

OUT_DIR = PROJECT / 'out'
OUT_DIR.mkdir(exist_ok=True, parents=True)

# =================== READ DATA ===================
df = pd.read_csv(CSV_PATH)
print(f"Read {len(df)} rows from {CSV_PATH.name}")

# =================== SUBCLASS ORDER (from CSV) ===================
subclass_meta = (
    df.drop_duplicates('subclass')
    .sort_values('subclass_order')[['subclass', 'broad_class', 'subclass_order']]
)
panel_subclasses = {
    bc: subclass_meta.loc[subclass_meta['broad_class'] == bc, 'subclass'].tolist()
    for bc in BROAD_CLASS_ORDER
}
subclass_order = [s for bc in BROAD_CLASS_ORDER for s in panel_subclasses[bc]]

# =================== GENE ORDER ===================
present_genes = set(df['gene'].unique())
gene_panel_map = {g: bc for g, bc in get_gene_panel_map().items() if g in present_genes}
gene_order = build_gene_order(gene_panel_map)
n_genes = len(gene_order)

print(f"  {len(subclass_order)} subclasses, {n_genes} genes")
for bc in BROAD_CLASS_ORDER:
    n_sub = len(panel_subclasses[bc])
    n_gen = sum(1 for g in gene_order if gene_panel_map.get(g) == bc)
    print(f"    {bc}: {n_sub} subclasses, {n_gen} genes")

# =================== PIVOT ===================
mean_piv = df.pivot(index='gene', columns='subclass', values='mean_expr').reindex(
    index=gene_order, columns=subclass_order
)
pct_piv = df.pivot(index='gene', columns='subclass', values='pct_expr').reindex(
    index=gene_order, columns=subclass_order
)

if STD_SCALE:
    row_min = mean_piv.min(axis=1)
    row_max = mean_piv.max(axis=1)
    denom = (row_max - row_min).replace(0, 1)
    mean_piv = mean_piv.sub(row_min, axis=0).div(denom, axis=0)

pct_capped = pct_piv.clip(upper=DOT_MAX)

# =================== COLOURS ===================
subclass_colours = load_subclass_colors(COLOR_CSV, subclasses=subclass_order)
vmin = float(np.nanmin(mean_piv.values))
vmax = float(np.nanmax(mean_piv.values))

# =================== FIGURE ===================
panel_widths = [max(len(panel_subclasses[bc]), 1) for bc in BROAD_CLASS_ORDER]
fig, axes = plt.subplots(
    1, 3, sharey=True,
    figsize=(FIG_WIDTH, FIG_HEIGHT),
    gridspec_kw={'width_ratios': panel_widths, 'wspace': 0.05},
)

sc_for_cbar = None
for panel_idx, bc in enumerate(BROAD_CLASS_ORDER):
    ax = axes[panel_idx]
    subs = panel_subclasses[bc]
    if not subs:
        ax.set_visible(False)
        continue

    sc = draw_dotplot_panel(
        ax, gene_order, subs, mean_piv, pct_capped,
        vmin, vmax, subclass_colours,
        expr_cmap=EXPR_CMAP, dot_scale=DOT_SCALE,
    )
    if sc_for_cbar is None:
        sc_for_cbar = sc

    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels(subs, rotation=90, ha='center', fontsize=7)
    ax.set_title(bc, fontsize=11, fontweight='bold', pad=8)

    if panel_idx == 0:
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(gene_order, fontsize=7)
    else:
        ax.tick_params(left=False)

cbar = fig.colorbar(
    sc_for_cbar, ax=axes.tolist(),
    shrink=0.25, pad=0.02, aspect=15, location='right',
)
cbar.set_label(
    'Mean expression\n(scaled per gene)' if STD_SCALE else 'Mean expression',
    fontsize=9,
)

add_size_legend(axes[-1], DOT_SCALE, DOT_MAX)
fig.suptitle('Marker gene dotplot', fontsize=13, fontweight='bold', y=0.98)

# =================== SAVE ===================
pdf_path = OUT_DIR / f'marker_dotplot_{stamp}.pdf'
png_path = OUT_DIR / f'marker_dotplot_{stamp}.png'
fig.savefig(pdf_path, bbox_inches='tight')
fig.savefig(png_path, dpi=300, bbox_inches='tight')
print(f"Saved {pdf_path.name}  and  {png_path.name}")
print(f"\n------ 02_plot_markers completed ------")
