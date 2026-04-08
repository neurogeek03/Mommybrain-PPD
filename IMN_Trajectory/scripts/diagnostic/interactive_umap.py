"""
Interactive UMAP with dropdown to switch between:
  - Cell type (subclass_name) — discrete colors + legend
  - Dcx, Ncam1, Rbfox3, Mki67 expression — continuous YlOrRd colorscale

Hover shows: cell type, sample, treatment.

Usage:
    conda run -n anndata_env python scripts/interactive_umap.py
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import plotly.graph_objects as go
import plotly.colors as pc

# ── paths ──────────────────────────────────────────────────────────────────────
DATA = "data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
OUT_DIR = "out"
os.makedirs(OUT_DIR, exist_ok=True)

# ── load ───────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(DATA)
print(f"  Shape: {adata.shape}")

# ── locate gene_symbol column ──────────────────────────────────────────────────
sym_cols = [c for c in adata.var.columns if "gene_symbol" in c.lower()]
if not sym_cols:
    raise ValueError(f"No gene_symbol column found. Available: {adata.var.columns.tolist()}")
sym_col = sym_cols[0]
print(f"  Using var column: '{sym_col}'")
symbols = adata.var[sym_col].astype(str)

# ── gene list (BrdU skipped; PSA-NCAM → Ncam1; NeuN → Rbfox3; Ki67 → Mki67) ──
TARGET_GENES = ["Dcx", "Ncam1", "Rbfox3", "Mki67"]

def extract_gene(adata, symbols, gene_name):
    """Return (display_label, expr_array) or (None, None) if not found."""
    mask = symbols.str.lower() == gene_name.lower()
    hits = adata.var_names[mask].tolist()
    if not hits:
        mask = symbols.str.lower().str.contains(gene_name.lower(), na=False)
        hits = adata.var_names[mask].tolist()
    if not hits:
        return None, None
    var_id = hits[0]
    label = symbols[var_id]
    idx = list(adata.var_names).index(var_id)
    expr = adata.X[:, idx]
    if sp.issparse(expr):
        expr = np.asarray(expr.todense()).flatten()
    else:
        expr = np.asarray(expr).flatten()
    return label, expr

# ── check UMAP ─────────────────────────────────────────────────────────────────
if "X_umap" not in adata.obsm:
    raise ValueError("X_umap not found in adata.obsm. Compute UMAP first.")

umap = adata.obsm["X_umap"]

# ── build base DataFrame ───────────────────────────────────────────────────────
df = pd.DataFrame({
    "UMAP1":        umap[:, 0],
    "UMAP2":        umap[:, 1],
    "subclass_name": adata.obs["subclass_name"].astype(str).values,
    "sample":        adata.obs["sample"].astype(str).values,
    "treatment":     adata.obs["treatment"].astype(str).values,
})

found_genes = {}  # original_name → display_label
for g in TARGET_GENES:
    label, expr = extract_gene(adata, symbols, g)
    if label is not None:
        df[label] = expr
        found_genes[g] = label
        print(f"  Found: {label}  (max={expr.max():.2f}, n_pos={( expr>0).sum()})")
    else:
        print(f"  NOT FOUND: {g}")

# ── categorical color map for cell types ───────────────────────────────────────
sorted_cats = sorted(df["subclass_name"].unique())
palette = (pc.qualitative.Plotly + pc.qualitative.D3 +
           pc.qualitative.G10 + pc.qualitative.Pastel)
cat_color = {cat: palette[i % len(palette)] for i, cat in enumerate(sorted_cats)}

HOVER_TMPL = (
    "<b>%{customdata[0]}</b><br>"
    "Sample: %{customdata[1]}<br>"
    "Treatment: %{customdata[2]}"
    "<extra></extra>"
)
MARKER_SIZE = 3

# ── build traces ───────────────────────────────────────────────────────────────
traces = []

# — cell-type traces (one per category for legend) —
for cat in sorted_cats:
    m = df["subclass_name"] == cat
    traces.append(go.Scatter(
        x=df.loc[m, "UMAP1"],
        y=df.loc[m, "UMAP2"],
        mode="markers",
        name=cat,
        marker=dict(color=cat_color[cat], size=MARKER_SIZE, opacity=0.85, line=dict(width=0)),
        customdata=df.loc[m, ["subclass_name", "sample", "treatment"]].values,
        hovertemplate=HOVER_TMPL,
        visible=True,
        legendgroup="celltype",
        showlegend=True,
    ))

n_cat_traces = len(sorted_cats)

# — gene expression traces (one per gene) —
for g_orig, g_label in found_genes.items():
    traces.append(go.Scatter(
        x=df["UMAP1"],
        y=df["UMAP2"],
        mode="markers",
        name=g_label,
        marker=dict(
            color=df[g_label],
            colorscale="YlOrRd",
            size=MARKER_SIZE,
            opacity=0.85,
            line=dict(width=0),
            showscale=True,
            colorbar=dict(
                title=dict(text=g_label, font=dict(color="black")),
                tickfont=dict(color="black"),
                outlinecolor="black",
                outlinewidth=0.5,
            ),
        ),
        customdata=df[["subclass_name", "sample", "treatment"]].values,
        hovertemplate=HOVER_TMPL,
        visible=False,
        showlegend=False,
    ))

n_gene_traces = len(found_genes)

# ── dropdown buttons ───────────────────────────────────────────────────────────
def make_vis(show_celltype=False, gene_idx=None):
    vis = [show_celltype] * n_cat_traces
    for i in range(n_gene_traces):
        vis.append(i == gene_idx if gene_idx is not None else False)
    return vis

buttons = [dict(
    label="Cell Type",
    method="update",
    args=[{"visible": make_vis(show_celltype=True)}, {"showlegend": True}],
)]

for i, (g_orig, g_label) in enumerate(found_genes.items()):
    buttons.append(dict(
        label=g_label,
        method="update",
        args=[{"visible": make_vis(gene_idx=i)}, {"showlegend": False}],
    ))

# ── assemble figure ────────────────────────────────────────────────────────────
fig = go.Figure(data=traces)

fig.update_layout(
    updatemenus=[dict(
        buttons=buttons,
        direction="down",
        showactive=True,
        x=1.01,
        xanchor="left",
        y=1.0,
        yanchor="top",
        bgcolor="white",
        bordercolor="black",
        borderwidth=1,
        font=dict(color="black"),
    )],
    plot_bgcolor="white",
    paper_bgcolor="white",
    font=dict(color="black", family="Arial", size=12),
    xaxis=dict(
        title="UMAP 1",
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="black",
        ticks="outside",
        tickcolor="black",
    ),
    yaxis=dict(
        title="UMAP 2",
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="black",
        ticks="outside",
        tickcolor="black",
    ),
    legend=dict(
        title="Cell Type",
        itemsizing="constant",
        font=dict(color="black"),
        bordercolor="black",
        borderwidth=0.5,
    ),
    title=dict(text="IMN Trajectory — UMAP", font=dict(color="black", size=14)),
    width=1000,
    height=750,
    margin=dict(r=180, t=60),
)

# ── save ───────────────────────────────────────────────────────────────────────
out_path = os.path.join(OUT_DIR, "interactive_umap.html")
fig.write_html(out_path, include_plotlyjs="cdn")
print(f"Saved: {out_path}")
