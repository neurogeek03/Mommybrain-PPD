"""
Marker gene visualization within '037 DG Glut' split by condition.

Stage 1 of DG trajectory analysis — gating check before pseudotime.
See docs/dg_trajectory_rationale.md for rationale and interpretation guide.

Produces:
  out/markers_037_umap_{stamp}.pdf/.png   — 6 genes x 2 conditions UMAP grid
  out/markers_037_violin_{stamp}.pdf/.png — violin per condition

Usage:
    conda run -n anndata_env python scripts/diagnostic/plot_markers_037.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from pathlib import Path
from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# ── config ───────────────────────────────────────────────────────────────────────
project_path    = Path(__file__).resolve().parents[2]
DATA            = project_path / 'data' / 'PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'
OUT_DIR         = project_path / 'out'
OUT_DIR.mkdir(exist_ok=True, parents=True)

SUBCLASS_TARGET = '037 DG Glut'
SUBCLASS_COL    = 'subclass_name'
SAMPLE_COL      = 'sample'
GENE_SYMBOL_COL = 'gene_symbol'

CORT_SAMPLES    = {'BC28', 'BC3', 'BC9'}
OIL_SAMPLES     = {'BC15', 'BC14', 'BC13'}

MARKERS         = ['Dcx', 'Sox2', 'Prox1', 'Calb1', 'Hopx', 'Mki67']

CONDITION_COLORS = {'CORT': '#2ca02c', 'OIL': '#ff7f0e'}

# ── load ─────────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(DATA)
print(f"  Full dataset: {adata.shape}")

# ── filter to 037 DG Glut ────────────────────────────────────────────────────────
if SUBCLASS_TARGET not in adata.obs[SUBCLASS_COL].unique():
    raise ValueError(f"'{SUBCLASS_TARGET}' not found in {SUBCLASS_COL}")
adata = adata[adata.obs[SUBCLASS_COL] == SUBCLASS_TARGET].copy()
print(f"  Filtered to '{SUBCLASS_TARGET}': {adata.n_obs} cells")

# ── condition label ───────────────────────────────────────────────────────────────
def assign_condition(sample):
    if sample in CORT_SAMPLES:
        return 'CORT'
    elif sample in OIL_SAMPLES:
        return 'OIL'
    return 'unknown'

adata.obs['condition'] = adata.obs[SAMPLE_COL].map(assign_condition)
print(f"  Condition counts:\n{adata.obs['condition'].value_counts().to_string()}")

# ── set var_names to gene symbols ────────────────────────────────────────────────
if GENE_SYMBOL_COL in adata.var.columns:
    adata.var_names = adata.var[GENE_SYMBOL_COL].astype(str)
    adata.var_names_make_unique()
    print(f"  var_names set to '{GENE_SYMBOL_COL}'")
else:
    print(f"  WARNING: '{GENE_SYMBOL_COL}' not in adata.var — keeping existing var_names")

# ── check which markers are present ──────────────────────────────────────────────
var_names_lower = {v.lower(): v for v in adata.var_names}
available = []
for gene in MARKERS:
    canonical = var_names_lower.get(gene.lower())
    if canonical:
        available.append(canonical)
        print(f"  Found: {canonical}")
    else:
        print(f"  WARNING: '{gene}' not found — skipping")
MARKERS = available

# ── normalize (log1p) for visualization ──────────────────────────────────────────
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
print("  Normalized + log1p.")

# ── extract expression for plotting before scaling clobbers X ────────────────────
def get_expr(gene):
    idx = list(adata.var_names).index(gene)
    col = adata.X[:, idx]
    return np.asarray(col.todense()).flatten() if sp.issparse(col) else np.asarray(col).flatten()

expr_cache = {gene: get_expr(gene) for gene in MARKERS}
print("  Expression cached pre-scaling.")

# ── HVGs → scale → PCA → neighbors → UMAP ───────────────────────────────────────
print("Computing UMAP on 037 cells...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
print("  UMAP done.")

umap       = adata.obsm['X_umap']
conditions = adata.obs['condition'].values
cond_order = ['CORT', 'OIL']

# ── PLOT 1: UMAP grid (n_genes × 2 conditions) ───────────────────────────────────
print("Plotting UMAP grid...")

n_genes = len(MARKERS)
fig, axes = plt.subplots(n_genes, 2, figsize=(10, n_genes * 3.5))
if n_genes == 1:
    axes = axes[np.newaxis, :]

for row, gene in enumerate(MARKERS):
    expr    = expr_cache[gene]
    pos_val = expr[expr > 0]
    vmax    = float(np.percentile(pos_val, 95)) if len(pos_val) > 0 else 1.0
    vmax    = max(vmax, 0.01)

    pct_pos = (expr > 0).mean() * 100
    print(f"  {gene}: {pct_pos:.1f}% cells expressing, vmax(p95)={vmax:.2f}")

    for col, cond in enumerate(cond_order):
        ax         = axes[row, col]
        mask_cond  = conditions == cond
        mask_other = ~mask_cond

        ax.scatter(umap[mask_other, 0], umap[mask_other, 1],
                   c='#e0e0e0', s=1, linewidths=0, rasterized=True, zorder=1)

        sc_plot = ax.scatter(umap[mask_cond, 0], umap[mask_cond, 1],
                             c=expr[mask_cond], cmap='YlOrRd',
                             vmin=0, vmax=vmax,
                             s=2, linewidths=0, rasterized=True, zorder=2)

        if col == 0:
            ax.set_ylabel(gene, fontsize=10, fontweight='bold')
        if row == 0:
            ax.set_title(cond, fontsize=11,
                         color=CONDITION_COLORS[cond], fontweight='bold')

        ax.axis('off')
        plt.colorbar(sc_plot, ax=ax, shrink=0.7, pad=0.02, label='log1p expr')

fig.suptitle(f"Marker genes in '037 DG Glut' — CORT vs OIL", fontsize=12, y=1.01)
fig.tight_layout()

umap_pdf = OUT_DIR / f"markers_037_umap_{stamp}.pdf"
umap_png = OUT_DIR / f"markers_037_umap_{stamp}.png"
fig.savefig(umap_pdf, bbox_inches='tight')
fig.savefig(umap_png, bbox_inches='tight', dpi=150)
plt.close(fig)
print(f"Saved: {umap_pdf.name}")

# ── PLOT 2: violin plots per condition ────────────────────────────────────────────
print("Plotting violins...")

fig2, axes2 = plt.subplots(1, n_genes, figsize=(n_genes * 3, 4))
if n_genes == 1:
    axes2 = [axes2]

for ax, gene in zip(axes2, MARKERS):
    expr         = expr_cache[gene]
    data_by_cond = [expr[conditions == c] for c in cond_order]

    parts = ax.violinplot(data_by_cond, positions=range(len(cond_order)),
                          showmedians=True, showextrema=False)

    for pc, cond in zip(parts['bodies'], cond_order):
        pc.set_facecolor(CONDITION_COLORS[cond])
        pc.set_alpha(0.7)
    parts['cmedians'].set_color('black')
    parts['cmedians'].set_linewidth(1.5)

    ax.set_xticks(range(len(cond_order)))
    ax.set_xticklabels(cond_order, fontsize=9)
    ax.set_title(gene, fontsize=10, fontweight='bold')
    if ax is axes2[0]:
        ax.set_ylabel('log1p expr', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig2.suptitle(f"Expression distribution — '037 DG Glut'", fontsize=11)
fig2.tight_layout()

vln_pdf = OUT_DIR / f"markers_037_violin_{stamp}.pdf"
vln_png = OUT_DIR / f"markers_037_violin_{stamp}.png"
fig2.savefig(vln_pdf, bbox_inches='tight')
fig2.savefig(vln_png, bbox_inches='tight', dpi=150)
plt.close(fig2)
print(f"Saved: {vln_pdf.name}")

# ── expression summary ────────────────────────────────────────────────────────────
print("\n── Expression summary (% cells > 0) ──")
for cond in cond_order:
    mask_c = conditions == cond
    n_cond = int(mask_c.sum())
    print(f"\n  {cond} (n={n_cond}):")
    for gene in MARKERS:
        expr     = expr_cache[gene]
        pct      = (expr[mask_c] > 0).mean() * 100
        pos_vals = expr[mask_c & (expr > 0)]
        mean_pos = float(pos_vals.mean()) if len(pos_vals) > 0 else 0.0
        print(f"    {gene:<10} {pct:5.1f}% positive,  mean in pos cells = {mean_pos:.3f}")

print(f"\n------ Script completed successfully ------")
