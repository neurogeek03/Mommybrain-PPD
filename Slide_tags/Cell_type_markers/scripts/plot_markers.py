"""
Panelled dotplot from marker_dotplot_data CSV.

Reads the CSV exported by create_marker_csv.py and produces a 3-panel dotplot
(Glutamatergic | GABAergic | Non-neuronal) with diagonal gene ordering per panel,
coloured subclass bar under the x-axis, and shared y-axis gene labels.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ plot_markers started at {stamp} ------")


# =================== PARAMS ===================
project_path = Path(__file__).resolve().parents[1]
CSV_PATH   = project_path / 'out' / 'marker_dotplot_data.csv'
COLOR_CSV  = Path('/scratch/mfafouti/Mommybrain/cluster_annotation_term.csv')

FIG_WIDTH   = 10
FIG_HEIGHT  = 14
EXPR_CMAP   = 'Purples'
DOT_MAX     = 0.8        # cap pct_expr for dot sizing (fraction)
DOT_SCALE   = 200        # base size multiplier for dots
STD_SCALE   = True       # standardise colour per gene (like scanpy standard_scale='var')

output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

BROAD_CLASS_ORDER = ['Glutamatergic', 'GABAergic', 'Non-neuronal']

# =================== READ DATA ===================

if not CSV_PATH:
    raise SystemExit("ERROR: Set CSV_PATH at the top of the script before running.")

df = pd.read_csv(CSV_PATH)
print(f"Read {len(df)} rows from {CSV_PATH}")

# =================== LOAD SUBCLASS COLOURS ===================

color_df = pd.read_csv(COLOR_CSV, usecols=['name', 'color_hex_triplet'])
color_df['name'] = color_df['name'].str.replace(r'[ /-]', '_', regex=True)
label_to_hex = dict(zip(color_df['name'], color_df['color_hex_triplet']))

def _normalise_name(s):
    """Replace spaces, hyphens, slashes with underscores for colour matching."""
    import re
    return re.sub(r'[ /-]', '_', s)

# =================== RECONSTRUCT ORDER ===================

df = df.sort_values(['subclass_order', 'gene'])
subclass_meta = df.drop_duplicates('subclass').sort_values('subclass_order')[
    ['subclass', 'broad_class', 'subclass_order']
]

# Build per-broad-class subclass lists (preserving existing sort)
panel_subclasses = {}
for bc in BROAD_CLASS_ORDER:
    panel_subclasses[bc] = subclass_meta.loc[
        subclass_meta['broad_class'] == bc, 'subclass'
    ].tolist()

subclass_order = []
for bc in BROAD_CLASS_ORDER:
    subclass_order.extend(panel_subclasses[bc])

# Gene owner = subclass where it has the highest specificity (i.e. the
# subclass it was selected as a marker for).
subclass_to_bc = dict(zip(subclass_meta['subclass'], subclass_meta['broad_class']))
subclass_to_order = dict(zip(subclass_meta['subclass'], subclass_meta['subclass_order']))

best_idx = df.groupby('gene')['specificity'].idxmax()
gene_owner_subclass = dict(zip(df.loc[best_idx, 'gene'], df.loc[best_idx, 'subclass']))
gene_owner_bc = {g: subclass_to_bc[s] for g, s in gene_owner_subclass.items()}

# Reorder genes: group by panel, then sort within panel by owner's
# subclass_order so that each panel's markers fall along its diagonal.
gene_order = []
for bc in BROAD_CLASS_ORDER:
    bc_genes = [g for g, obc in gene_owner_bc.items() if obc == bc]
    bc_genes.sort(key=lambda g: subclass_to_order[gene_owner_subclass[g]])
    gene_order.extend(bc_genes)

print(f"  {len(subclass_order)} subclasses, {len(gene_order)} genes")
for bc in BROAD_CLASS_ORDER:
    n_sub = len(panel_subclasses[bc])
    n_gen = sum(1 for g in gene_order if gene_owner_bc[g] == bc)
    print(f"    {bc}: {n_sub} subclasses, {n_gen} genes")

# =================== PIVOT ===================

mean_piv = df.pivot(index='gene', columns='subclass', values='mean_expr').reindex(
    index=gene_order, columns=subclass_order
)
pct_piv = df.pivot(index='gene', columns='subclass', values='pct_expr').reindex(
    index=gene_order, columns=subclass_order
)

# Optional per-gene standardisation (0-1 scale per row)
if STD_SCALE:
    row_min = mean_piv.min(axis=1)
    row_max = mean_piv.max(axis=1)
    denom = (row_max - row_min).replace(0, 1)
    mean_piv = mean_piv.sub(row_min, axis=0).div(denom, axis=0)

# Cap pct for sizing
pct_capped = pct_piv.clip(upper=DOT_MAX)

# =================== COLOUR MAP FOR SUBCLASSES ===================

subclass_colours = {}
for sub in subclass_order:
    key = _normalise_name(sub)
    subclass_colours[sub] = label_to_hex.get(key, '#808080')

# =================== PLOT (3 PANELS) ===================

n_genes = len(gene_order)
panel_widths = [len(panel_subclasses[bc]) for bc in BROAD_CLASS_ORDER]

fig, axes = plt.subplots(
    1, 3, sharey=True,
    figsize=(FIG_WIDTH, FIG_HEIGHT),
    gridspec_kw={'width_ratios': panel_widths, 'wspace': 0.05}
)

# We need a consistent colour scale across panels
vmin = mean_piv.values[~np.isnan(mean_piv.values)].min()
vmax = mean_piv.values[~np.isnan(mean_piv.values)].max()

BAR_HEIGHT = 0.8  # height of the colour bar rectangles
bar_y = n_genes + 0.5  # y position for colour bar (below gene grid, since y is inverted)

sc_for_cbar = None  # will store one scatter for the colourbar

for panel_idx, bc in enumerate(BROAD_CLASS_ORDER):
    ax = axes[panel_idx]
    subs = panel_subclasses[bc]
    n_subs = len(subs)

    # Scatter: all genes (y) vs this panel's subclasses (x)
    xs, ys, sizes, colours = [], [], [], []
    for gi, gene in enumerate(gene_order):
        for si, sub in enumerate(subs):
            xs.append(si)
            ys.append(gi)
            sizes.append(pct_capped.at[gene, sub] * DOT_SCALE)
            colours.append(mean_piv.at[gene, sub])

    sc = ax.scatter(xs, ys, s=sizes, c=colours, cmap=EXPR_CMAP,
                    vmin=vmin, vmax=vmax, edgecolors='none')
    if sc_for_cbar is None:
        sc_for_cbar = sc

    # --- Coloured x-axis bar ---
    for si, sub in enumerate(subs):
        hex_col = subclass_colours[sub]
        rect = Rectangle(
            (si - 0.4, bar_y), 0.8, BAR_HEIGHT,
            facecolor=hex_col, edgecolor='none', clip_on=False
        )
        ax.add_patch(rect)

    # --- Axis formatting ---
    ax.set_xlim(-0.5, n_subs - 0.5)
    ax.set_ylim(bar_y + BAR_HEIGHT + 0.2, -0.5)  # inverted y; extend to show bar

    ax.set_xticks(range(n_subs))
    ax.set_xticklabels(subs, rotation=90, ha='center', fontsize=7)

    if panel_idx == 0:
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(gene_order, fontsize=7)
    else:
        ax.tick_params(left=False)

    ax.set_title(bc, fontsize=11, fontweight='bold', pad=8)

    # Light grid lines for readability
    for gi in range(n_genes):
        ax.axhline(gi, color='#f0f0f0', linewidth=0.3, zorder=0)

# --- Colourbar (mean expression) ---
cbar = fig.colorbar(sc_for_cbar, ax=axes.tolist(), shrink=0.25, pad=0.02,
                    aspect=15, location='right')
cbar.set_label('Mean expression\n(scaled per gene)' if STD_SCALE else 'Mean expression',
               fontsize=9)

# --- Size legend (pct expressing) ---
legend_pcts = [0.2, 0.4, 0.6, 0.8]
legend_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
           markersize=np.sqrt(p * DOT_SCALE) * 0.6, label=f'{int(p*100)}%')
    for p in legend_pcts if p <= DOT_MAX
]
axes[-1].legend(handles=legend_handles, title='% expressing', loc='upper left',
                bbox_to_anchor=(1.12, 0.85), frameon=False, fontsize=8,
                title_fontsize=9)

fig.suptitle('Marker gene dotplot', fontsize=13, fontweight='bold', y=0.98)

# =================== SAVE ===================

pdf_path = output_base / f"marker_dotplot_{stamp}.pdf"
png_path = output_base / f"marker_dotplot_{stamp}.png"
fig.savefig(pdf_path, bbox_inches='tight')
fig.savefig(png_path, dpi=300, bbox_inches='tight')
print(f"Saved {pdf_path.name} and {png_path.name}")

print(f"\n------ plot_markers completed successfully ------")
