"""
PCA on HVGs for 3 IMN-related subclasses:
  - 037 DG Glut
  - 038 DG-PIR Ex IMN
  - 045 OB-STR-CTX Inh IMN

Produces:
  out/pca_hvg_dg_imn_variance_{stamp}.pdf   — scree plot (variance explained per PC)
  out/pca_hvg_dg_imn_loadings_{stamp}.pdf   — top gene loadings per PC (heatmap)
  out/pca_hvg_dg_imn_loadings_{stamp}.csv   — loadings table

Usage:
    conda run -n anndata_env python scripts/pca_hvg_dg_imn.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# ── config ──────────────────────────────────────────────────────────────────────
project_path = Path.cwd().parents[0]
DATA         = project_path / 'data' / 'PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'
OUT_DIR      = project_path / 'out'
OUT_DIR.mkdir(exist_ok=True, parents=True)

SUBCLASS_COL    = "subclass_name"
SAMPLE_COL      = "sample"
GENE_SYMBOL_COL = "gene_symbol"
N_HVG           = 2000
N_PCS           = 50
TOP_N_GENES     = 15   # top loading genes to show per PC in the heatmap

TARGET_TYPES = [
    '037 DG Glut',
    '038 DG-PIR Ex IMN',
    '045 OB-STR-CTX Inh IMN',
]

COLORS = {
    '037 DG Glut':           '#4C72B0',
    '038 DG-PIR Ex IMN':     '#DD8452',
    '045 OB-STR-CTX Inh IMN': '#55A868',
}

# ── load ────────────────────────────────────────────────────────────────────────
adata = sc.read_h5ad(DATA)
print(adata)

# ── filter to target subclasses ─────────────────────────────────────────────────
missing = [ct for ct in TARGET_TYPES if ct not in adata.obs[SUBCLASS_COL].unique()]
if missing:
    raise ValueError(f"Subclasses not found: {missing}")

adata = adata[adata.obs[SUBCLASS_COL].isin(TARGET_TYPES)].copy()
adata.obs[SUBCLASS_COL] = pd.Categorical(adata.obs[SUBCLASS_COL], categories=TARGET_TYPES)
print(f"Filtered to {adata.n_obs} cells")
print(adata.obs[SUBCLASS_COL].value_counts().to_string())

# ── cell counts per subclass × sample ───────────────────────────────────────────
counts = (
    adata.obs.groupby([SUBCLASS_COL, SAMPLE_COL], observed=True)
    .size()
    .unstack(fill_value=0)
)
counts['TOTAL'] = counts.sum(axis=1)
print(f"\nCell counts per subclass × sample:")
print(counts.to_string())

counts_csv = OUT_DIR / f"pca_hvg_dg_imn_cell_counts_{stamp}.csv"
counts.to_csv(counts_csv)
print(f"Saved: {counts_csv.name}")

# ── gene symbols ────────────────────────────────────────────────────────────────
if GENE_SYMBOL_COL in adata.var.columns:
    adata.var_names = adata.var[GENE_SYMBOL_COL].astype(str)
    adata.var_names_make_unique()
    print(f"var_names set to '{GENE_SYMBOL_COL}'")
else:
    print(f"WARNING: '{GENE_SYMBOL_COL}' not in adata.var — keeping existing var_names")

# # ── normalize ───────────────────────────────────────────────────────────────────
# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata)
# print("Normalize + log1p done.")

# ── HVGs ────────────────────────────────────────────────────────────────────────
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, subset=True)
print(f"HVGs selected: {adata.n_vars}")

# ── scale & PCA ─────────────────────────────────────────────────────────────────
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=N_PCS, use_highly_variable=False)
print("PCA done.")

var_ratio = adata.uns['pca']['variance_ratio']   # shape (N_PCS,)
loadings  = adata.varm['PCs']                    # shape (n_hvg, N_PCS)

# ── PLOT 1: scree plot ────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# individual variance explained
axes[0].bar(range(1, N_PCS + 1), var_ratio * 100, color='steelblue', edgecolor='none')
axes[0].set_xlabel('PC')
axes[0].set_ylabel('Variance explained (%)')
axes[0].set_title('Variance explained per PC')
axes[0].set_xlim(0.5, N_PCS + 0.5)

# cumulative
axes[1].plot(range(1, N_PCS + 1), np.cumsum(var_ratio) * 100, marker='o', markersize=3,
             color='steelblue')
axes[1].axhline(80, color='red', linestyle='--', linewidth=0.8, label='80%')
axes[1].axhline(90, color='orange', linestyle='--', linewidth=0.8, label='90%')
axes[1].set_xlabel('PC')
axes[1].set_ylabel('Cumulative variance explained (%)')
axes[1].set_title('Cumulative variance explained')
axes[1].set_xlim(0.5, N_PCS + 0.5)
axes[1].legend(fontsize=8)

fig.suptitle('PCA on HVGs — DG Glut / DG-PIR Ex IMN / OB-STR-CTX Inh IMN', fontsize=11)
fig.tight_layout()

var_path_pdf = OUT_DIR / f"pca_hvg_dg_imn_variance_{stamp}.pdf"
var_path_png = OUT_DIR / f"pca_hvg_dg_imn_variance_{stamp}.png"
fig.savefig(var_path_pdf, bbox_inches='tight')
fig.savefig(var_path_png, bbox_inches='tight', dpi=150)
plt.close(fig)
print(f"Saved: {var_path_pdf.name}")

# ── PLOT 2: top gene loadings per PC (heatmap) ───────────────────────────────
# Show TOP_N_GENES with largest absolute loading for each of the first N_PCS_SHOW PCs
N_PCS_SHOW = min(20, N_PCS)

gene_names = adata.var_names.tolist()

# collect top genes per PC and their loadings
records = []
for pc_idx in range(N_PCS_SHOW):
    pc_loadings = loadings[:, pc_idx]
    top_idx = np.argsort(np.abs(pc_loadings))[::-1][:TOP_N_GENES]
    for rank, g_idx in enumerate(top_idx):
        records.append({
            'PC':     pc_idx + 1,
            'gene':   gene_names[g_idx],
            'loading': pc_loadings[g_idx],
            'abs_loading': abs(pc_loadings[g_idx]),
            'rank':   rank + 1,
        })

loadings_df = pd.DataFrame(records)

# Build a union of all top genes across the shown PCs, then make a heatmap
top_genes_union = loadings_df.groupby('gene')['abs_loading'].max().nlargest(TOP_N_GENES * 2).index.tolist()

heatmap_data = pd.DataFrame(
    index=top_genes_union,
    columns=[f"PC{i+1}" for i in range(N_PCS_SHOW)],
    dtype=float,
)
for pc_idx in range(N_PCS_SHOW):
    pc_col = f"PC{pc_idx+1}"
    pc_loadings = loadings[:, pc_idx]
    for gene in top_genes_union:
        if gene in gene_names:
            g_idx = gene_names.index(gene)
            heatmap_data.loc[gene, pc_col] = pc_loadings[g_idx]

heatmap_data = heatmap_data.fillna(0).astype(float)
# sort genes by PC1 loading
heatmap_data = heatmap_data.sort_values('PC1', ascending=False)

fig2, ax = plt.subplots(figsize=(N_PCS_SHOW * 0.6 + 2, len(top_genes_union) * 0.4 + 1.5))
im = ax.imshow(heatmap_data.values, aspect='auto', cmap='RdBu_r',
               vmin=-heatmap_data.values.__abs__().max(),
               vmax= heatmap_data.values.__abs__().max())
ax.set_xticks(range(N_PCS_SHOW))
ax.set_xticklabels([f"PC{i+1}" for i in range(N_PCS_SHOW)], fontsize=7, rotation=45, ha='right')
ax.set_yticks(range(len(top_genes_union)))
ax.set_yticklabels(heatmap_data.index.tolist(), fontsize=7)
ax.set_title(f'Top gene loadings (PCs 1–{N_PCS_SHOW})', fontsize=10)
plt.colorbar(im, ax=ax, label='Loading', shrink=0.6)
fig2.tight_layout()

load_path_pdf = OUT_DIR / f"pca_hvg_dg_imn_loadings_{stamp}.pdf"
load_path_png = OUT_DIR / f"pca_hvg_dg_imn_loadings_{stamp}.png"
fig2.savefig(load_path_pdf, bbox_inches='tight')
fig2.savefig(load_path_png, bbox_inches='tight', dpi=150)
plt.close(fig2)
print(f"Saved: {load_path_pdf.name}")

# ── PLOT 3: bar charts of top +/- loading genes for first 10 PCs ─────────────
N_PCS_BAR = min(10, N_PCS)
fig3, axes3 = plt.subplots(2, 5, figsize=(20, 8)) if N_PCS_BAR == 10 else \
              plt.subplots(1, N_PCS_BAR, figsize=(4 * N_PCS_BAR, 4))
axes3 = np.array(axes3).flatten()

for pc_idx in range(N_PCS_BAR):
    ax = axes3[pc_idx]
    pc_loadings = loadings[:, pc_idx]
    top_pos_idx = np.argsort(pc_loadings)[::-1][:TOP_N_GENES]
    top_neg_idx = np.argsort(pc_loadings)[:TOP_N_GENES]
    combined_idx = np.concatenate([top_pos_idx, top_neg_idx])
    combined_idx = pd.unique(combined_idx)  # deduplicate

    genes_show = [gene_names[i] for i in combined_idx]
    vals_show  = [pc_loadings[i] for i in combined_idx]
    order = np.argsort(vals_show)
    genes_show = [genes_show[i] for i in order]
    vals_show  = [vals_show[i]  for i in order]

    colors = ['#d62728' if v > 0 else '#1f77b4' for v in vals_show]
    ax.barh(range(len(genes_show)), vals_show, color=colors)
    ax.set_yticks(range(len(genes_show)))
    ax.set_yticklabels(genes_show, fontsize=6)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_title(f'PC{pc_idx+1} ({var_ratio[pc_idx]*100:.1f}%)', fontsize=8)
    ax.set_xlabel('Loading', fontsize=7)

fig3.suptitle(f'Top {TOP_N_GENES} positive/negative loading genes per PC', fontsize=11)
fig3.tight_layout()

bar_path_pdf = OUT_DIR / f"pca_hvg_dg_imn_bar_{stamp}.pdf"
bar_path_png = OUT_DIR / f"pca_hvg_dg_imn_bar_{stamp}.png"
fig3.savefig(bar_path_pdf, bbox_inches='tight')
fig3.savefig(bar_path_png, bbox_inches='tight', dpi=150)
plt.close(fig3)
print(f"Saved: {bar_path_pdf.name}")

# ── save loadings CSV ────────────────────────────────────────────────────────
csv_path = OUT_DIR / f"pca_hvg_dg_imn_loadings_{stamp}.csv"
loadings_df.to_csv(csv_path, index=False)
print(f"Saved: {csv_path.name}")

# ── PLOT 4: PC6 scores on UMAP ───────────────────────────────────────────────
PC_TO_PLOT = 6  # 1-based

print(f"\nComputing neighbors + UMAP for PC{PC_TO_PLOT} overlay...")
sc.pp.neighbors(adata, n_pcs=N_PCS)
sc.tl.umap(adata)

pc_scores = adata.obsm['X_pca'][:, PC_TO_PLOT - 1]   # 0-based index
umap_coords = adata.obsm['X_umap']
subclass_labels = adata.obs[SUBCLASS_COL].astype(str).values

fig4, axes4 = plt.subplots(1, 2, figsize=(14, 5))

# panel A: PC6 scores as continuous colour
vmax = np.abs(pc_scores).max()
sc4 = axes4[0].scatter(
    umap_coords[:, 0], umap_coords[:, 1],
    c=pc_scores, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
    s=2, linewidths=0, rasterized=True,
)
plt.colorbar(sc4, ax=axes4[0], label=f'PC{PC_TO_PLOT} score', shrink=0.7)
axes4[0].set_title(f'PC{PC_TO_PLOT} scores on UMAP', fontsize=10)
axes4[0].axis('off')

# panel B: subclass labels for reference
for label, color in COLORS.items():
    m = subclass_labels == label
    axes4[1].scatter(
        umap_coords[m, 0], umap_coords[m, 1],
        c=color, s=2, linewidths=0, rasterized=True, label=label,
    )
axes4[1].set_title('Subclass', fontsize=10)
axes4[1].axis('off')
axes4[1].legend(markerscale=4, fontsize=7, loc='best', frameon=False)

fig4.suptitle(
    f'PC{PC_TO_PLOT} ({var_ratio[PC_TO_PLOT - 1] * 100:.1f}% variance) — DG Glut / DG-PIR Ex IMN / OB-STR-CTX Inh IMN',
    fontsize=10,
)
fig4.tight_layout()

umap_pdf = OUT_DIR / f"pca_hvg_dg_imn_pc{PC_TO_PLOT}_umap_{stamp}.pdf"
umap_png = OUT_DIR / f"pca_hvg_dg_imn_pc{PC_TO_PLOT}_umap_{stamp}.png"
fig4.savefig(umap_pdf, bbox_inches='tight')
fig4.savefig(umap_png, bbox_inches='tight', dpi=150)
plt.close(fig4)
print(f"Saved: {umap_pdf.name}")

print(f"\n------ Script completed successfully ------")