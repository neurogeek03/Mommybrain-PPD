import matplotlib
matplotlib.use('Agg')

from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path(__file__).resolve().parents[1]
output_base  = project_path / 'out' / 'IMN_clusters'
output_base.mkdir(exist_ok=True, parents=True)

# =================== PARAMS ===================
group_col       = "subclass_name"
GENE_SYMBOL_COL = 'gene_symbol'

TARGET = ['038 DG-PIR Ex IMN', '045 OB-STR-CTX Inh IMN']

# Curated markers for score-based relabeling
MARKERS_DG_PIR     = ['Prox1', 'Tbr1', 'Sox4', 'Pou3f3']              # 038 DG-PIR Ex IMN (excitatory TFs)
MARKERS_OB_STR_CTX = ['Pax6', 'Dlx2', 'Zeb2', 'Bcl11b', 'Meis2', 'Nfix']  # 045 OB-STR-CTX Inh IMN (inhibitory TFs)

# Additional markers shown in dotplot only (not used for scoring)
EXTRA_DG_PIR_MARKERS = ['Mex3a', 'Neurod1']  # 038 expression markers
DCX_GENE           = 'Dcx'                     # broad neurogenesis marker

# Broad excitatory/inhibitory markers for sanity check
EXCITATORY_MARKERS = ['Slc17a7', 'Slc17a6']
INHIBITORY_MARKERS = ['Gad1', 'Gad2']

# =================== INPUT ===================
ad_file = '/scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'

adata = sc.read_h5ad(ad_file)
print(adata)

# =================== GENE SYMBOLS ===================
if GENE_SYMBOL_COL in adata.var.columns:
    adata.var_names = adata.var[GENE_SYMBOL_COL].astype(str)
    adata.var_names_make_unique()
    print(f"\nvar_names set to: {GENE_SYMBOL_COL}")
else:
    print(f"\nWARNING: '{GENE_SYMBOL_COL}' not found — keeping existing var_names")
    print(f"  Available columns: {adata.var.columns.tolist()}")

# =================== SUBSET ===================
obs_unique = sorted(adata.obs[group_col].unique())
missing = [ct for ct in TARGET if ct not in obs_unique]
if missing:
    raise ValueError(f"Target cell types not found in '{group_col}': {missing}")

sub = adata[adata.obs[group_col].isin(TARGET)].copy()
print(f"\nSubset to {sub.n_obs} cells")
print(sub.obs[group_col].value_counts().to_string())

# # =================== NORMALIZE ===================
# sc.pp.normalize_total(sub)
# sc.pp.log1p(sub)
# print("\nNormalize + log1p done.")

# =================== CASE-INSENSITIVE GENE LOOKUP ===================
def resolve_gene(gene, var_names):
    """Return the exact var_name matching gene (case-insensitive). None if not found."""
    if gene in var_names:
        return gene
    matches = [g for g in var_names if g.lower() == gene.lower()]
    if matches:
        print(f"  Gene '{gene}' resolved to '{matches[0]}' (case-insensitive match)")
        return matches[0]
    print(f"  WARNING: '{gene}' not found in var_names — skipping")
    return None

def resolve_gene_list(gene_list, var_names, label=''):
    resolved = [resolve_gene(g, var_names) for g in gene_list]
    resolved = [g for g in resolved if g is not None]
    if label:
        print(f"  {label}: {resolved}")
    return resolved

print("\n--- Resolving gene names ---")
markers_dg_pir     = resolve_gene_list(MARKERS_DG_PIR,       sub.var_names, '038 DG-PIR Ex IMN (TFs)')
markers_ob_str_ctx = resolve_gene_list(MARKERS_OB_STR_CTX,  sub.var_names, '045 OB-STR-CTX Inh IMN (TFs)')
extra_dg_pir       = resolve_gene_list(EXTRA_DG_PIR_MARKERS, sub.var_names, '038 extra markers')
exc_markers        = resolve_gene_list(EXCITATORY_MARKERS,  sub.var_names, 'Excitatory (broad)')
inh_markers        = resolve_gene_list(INHIBITORY_MARKERS,  sub.var_names, 'Inhibitory (broad)')
dcx_resolved       = resolve_gene(DCX_GENE, sub.var_names)
dcx_list           = [dcx_resolved] if dcx_resolved else []

if not markers_dg_pir:
    raise RuntimeError("No markers resolved for 038 DG-PIR Ex IMN — cannot score.")
if not markers_ob_str_ctx:
    raise RuntimeError("No markers resolved for 045 OB-STR-CTX Inh IMN — cannot score.")

# =================== DIMENSIONALITY REDUCTION ===================
print("\nRunning HVG selection, PCA, neighbors, UMAP on subset...")
sc.pp.highly_variable_genes(sub, n_top_genes=2000)
sc.tl.pca(sub, use_highly_variable=True)
sc.pp.neighbors(sub, n_pcs=20)
sc.tl.umap(sub)

for res in [0.1, 0.3, 0.5]:
    sc.tl.leiden(sub, resolution=res, key_added=f'leiden_{res}')
    print(f"  leiden_{res}: {sub.obs[f'leiden_{res}'].nunique()} clusters")

# =================== UMAP — CURRENT LABELS + MARKERS ===================
import matplotlib.pyplot as plt

print("\nSaving UMAP visualizations...")

# Current labels
fig, ax = plt.subplots(figsize=(6, 5))
sc.pl.umap(sub, color='subclass_name', ax=ax, show=False)
fig.savefig(str(output_base / f"DG_imn_relabel_umap_labels_{stamp}.png"), bbox_inches='tight', dpi=150)
plt.close(fig)

# Leiden at each resolution
for res in [0.1, 0.3, 0.5]:
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(sub, color=f'leiden_{res}', ax=ax, show=False, title=f'Leiden {res}')
    fig.savefig(str(output_base / f"DG_imn_relabel_umap_leiden{res}_{stamp}.png"), bbox_inches='tight', dpi=150)
    plt.close(fig)

# Curated cluster markers
curated_umap_genes = markers_dg_pir + markers_ob_str_ctx + dcx_list
if curated_umap_genes:
    sc.pl.umap(sub, color=curated_umap_genes, ncols=2, show=False,
               save=f"_DG_imn_relabel_curated_markers_{stamp}.png")

# Broad excitatory/inhibitory markers
broad_umap_genes = exc_markers + inh_markers
if broad_umap_genes:
    sc.pl.umap(sub, color=broad_umap_genes, ncols=2, show=False,
               save=f"_DG_imn_relabel_broad_markers_{stamp}.png")

# Violin — broad markers by current label
if broad_umap_genes:
    sc.pl.violin(sub, keys=broad_umap_genes, groupby='subclass_name', show=False,
                 save=f"_DG_imn_relabel_violin_broad_{stamp}.png")

# =================== INTERACTIVE UMAP (Plotly) ===================
import plotly.graph_objects as go

# All marker genes to expose in the dropdown
viz_genes = list(dict.fromkeys(
    markers_dg_pir + extra_dg_pir + markers_ob_str_ctx + dcx_list + exc_markers + inh_markers
))

umap_x = sub.obsm['X_umap'][:, 0]
umap_y = sub.obsm['X_umap'][:, 1]
cell_labels = sub.obs['subclass_name'].astype(str).values

# Extract expression per gene (sparse-safe, no full matrix load)
expr_data = {}
for gene in viz_genes:
    col = sub[:, gene].X
    if sp.issparse(col):
        col = col.toarray().flatten()
    else:
        col = col.flatten()
    expr_data[gene] = col.astype(float)

# One scatter trace per gene; only the first is visible on load
fig_int = go.Figure()
for i, gene in enumerate(viz_genes):
    expr = expr_data[gene]
    fig_int.add_trace(go.Scatter(
        x=umap_x,
        y=umap_y,
        mode='markers',
        marker=dict(
            color=expr,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title=gene),
            size=5,
            opacity=0.8,
        ),
        text=[f"{c}<br>{gene}: {e:.3f}" for c, e in zip(cell_labels, expr)],
        hoverinfo='text',
        name=gene,
        visible=(i == 0),
    ))

# Dropdown: each button makes exactly one trace visible
buttons = [
    dict(
        label=gene,
        method='update',
        args=[
            {'visible': [j == i for j in range(len(viz_genes))]},
            {'title': f'UMAP — {gene}'},
        ],
    )
    for i, gene in enumerate(viz_genes)
]

fig_int.update_layout(
    title=f'UMAP — {viz_genes[0]}',
    xaxis_title='UMAP 1',
    yaxis_title='UMAP 2',
    updatemenus=[dict(
        active=0,
        buttons=buttons,
        direction='down',
        x=0.01,
        xanchor='left',
        y=1.12,
        yanchor='top',
        showactive=True,
    )],
    width=700,
    height=600,
)

html_path = output_base / f"DG_imn_relabel_interactive_{stamp}.html"
fig_int.write_html(str(html_path))
print(f"Saved: {html_path.name}")

# =================== DOTPLOT HELPER ===================
# Dotplot — curated + broad markers, current labels
dotplot_genes = {}
if markers_dg_pir:
    dotplot_genes['DG-PIR Ex IMN (TFs)'] = markers_dg_pir
if extra_dg_pir:
    dotplot_genes['DG-PIR Ex IMN (markers)'] = extra_dg_pir
if markers_ob_str_ctx:
    dotplot_genes['OB-STR-CTX Inh IMN (TFs)'] = markers_ob_str_ctx
if dcx_list:
    dotplot_genes['Neurogenesis'] = dcx_list
if exc_markers:
    dotplot_genes['Excitatory (broad)'] = exc_markers
if inh_markers:
    dotplot_genes['Inhibitory (broad)'] = inh_markers

dp = sc.pl.dotplot(
    sub,
    var_names=dotplot_genes,
    groupby='subclass_name',
    use_raw=False,
    standard_scale='var',
    swap_axes=True,
    return_fig=True,
    show=False,
)
main_ax = dp.get_axes()['mainplot_ax']
main_ax.set_xticklabels(main_ax.get_xticklabels(), rotation=45, ha='right')
fig = main_ax.get_figure()
fig.savefig(str(output_base / f"DG_imn_relabel_dotplot_before_{stamp}.pdf"), bbox_inches='tight')
fig.savefig(str(output_base / f"DG_imn_relabel_dotplot_before_{stamp}.png"), bbox_inches='tight', dpi=150)
plt.close(fig)
print(f"Saved: DG_imn_relabel_dotplot_before_{stamp}.pdf/png")

# =================== SCORE-BASED RELABELING ===================
print("\nComputing per-cell gene scores...")
sc.tl.score_genes(sub, gene_list=markers_dg_pir,     score_name='score_DG_PIR')
sc.tl.score_genes(sub, gene_list=markers_ob_str_ctx, score_name='score_OB_STR_CTX')

sub.obs['new_label'] = sub.obs.apply(
    lambda r: '038 DG-PIR Ex IMN' if r['score_DG_PIR'] >= r['score_OB_STR_CTX']
              else '045 OB-STR-CTX Inh IMN',
    axis=1
)

print("\nLabel counts BEFORE relabeling:")
print(sub.obs['subclass_name'].value_counts().to_string())
print("\nLabel counts AFTER relabeling:")
print(sub.obs['new_label'].value_counts().to_string())

# How many cells changed?
n_changed = (sub.obs['subclass_name'] != sub.obs['new_label']).sum()
pct_changed = n_changed / sub.n_obs * 100
print(f"\nCells relabeled: {n_changed} / {sub.n_obs} ({pct_changed:.1f}%)")

# UMAP — scores and new labels
sc.pl.umap(sub, color=['score_DG_PIR', 'score_OB_STR_CTX', 'new_label'], ncols=3, show=False,
           save=f"_DG_imn_relabel_scores_{stamp}.png")

# =================== VALIDATE — DOTPLOT WITH NEW LABELS ===================
dp2 = sc.pl.dotplot(
    sub,
    var_names=dotplot_genes,
    groupby='new_label',
    use_raw=False,
    standard_scale='var',
    swap_axes=True,
    return_fig=True,
    show=False,
)
main_ax2 = dp2.get_axes()['mainplot_ax']
main_ax2.set_xticklabels(main_ax2.get_xticklabels(), rotation=45, ha='right')
fig2 = main_ax2.get_figure()
fig2.savefig(str(output_base / f"DG_imn_relabel_dotplot_after_{stamp}.pdf"), bbox_inches='tight')
fig2.savefig(str(output_base / f"DG_imn_relabel_dotplot_after_{stamp}.png"), bbox_inches='tight', dpi=150)
plt.close(fig2)
print(f"Saved: DG_imn_relabel_dotplot_after_{stamp}.pdf/png")

# # =================== PROPAGATE LABELS TO FULL ADATA ===================
# print("\nPropagating new labels to full adata...")
# adata.obs.loc[sub.obs_names, group_col] = sub.obs['new_label']

# print("\nFull adata — counts for target subclasses after relabeling:")
# print(adata.obs.loc[adata.obs[group_col].isin(TARGET), group_col].value_counts().to_string())

# # =================== SAVE CORRECTED H5AD ===================
# out_h5ad = output_base / f'adata_relabeled_DGimn_{stamp}.h5ad'
# adata.write_h5ad(out_h5ad)
# print(f"\nSaved: {out_h5ad.name}")

# print(f"\n------ Script completed successfully ------")
