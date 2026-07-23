"""
Interactive UMAP for oligodendrocytes coloured by CORT-response gene expression
and UCell score. Outputs a self-contained HTML file.

Run with:
    uv run python code/03_oligo_umap.py
"""

import os
os.environ["NUMBA_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np
import anndata as ad
import scanpy as sc
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pyucell import compute_ucell_scores

H5AD = "data/slide_tags_129493.h5ad"
OUT_HTML = "data/oligo_umap.html"
CELL_TYPE_COL = "subclass_name"
TREATMENT_COL = "treatment"

CORT_GENES = [
    "Fkbp5", "Sgk1", "Tsc22d3", "Klf9", "Zbtb16", "Rcan2",
    "Hif3a", "Ddit4", "Errfi1", "Nfkbia", "Bcl6",
]

SIGNATURE = {"CORT_response": CORT_GENES}

MARKER_SIZE = 2
MARKER_OPACITY = 0.6


def log1p_expr(adata, gene):
    """Return log1p-normalised expression for a gene as a numpy array."""
    idx = adata.var_names.get_loc(gene)
    col = adata.X[:, idx]
    if hasattr(col, "toarray"):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return np.log1p(col)


def make_colorscale_trace(x, y, values, title, colorscale="Viridis", zmin=None, zmax=None):
    return go.Scattergl(
        x=x, y=y,
        mode="markers",
        marker=dict(
            size=MARKER_SIZE,
            opacity=MARKER_OPACITY,
            color=values,
            colorscale=colorscale,
            cmin=zmin if zmin is not None else float(np.nanmin(values)),
            cmax=zmax if zmax is not None else float(np.nanmax(values)),
            colorbar=dict(thickness=10, len=0.4),
            showscale=True,
        ),
        text=title,
        hovertemplate=f"<b>{title}</b><br>value: %{{marker.color:.3f}}<extra></extra>",
        name=title,
    )


def make_categorical_trace(x, y, labels, label_colors, title):
    traces = []
    for label, color in label_colors.items():
        mask = labels == label
        traces.append(go.Scattergl(
            x=x[mask], y=y[mask],
            mode="markers",
            marker=dict(size=MARKER_SIZE, opacity=MARKER_OPACITY, color=color),
            name=label,
            legendgroup=title,
            hovertemplate=f"<b>{label}</b><extra></extra>",
        ))
    return traces


def main():
    print("Loading AnnData...")
    adata = ad.read_h5ad(H5AD)
    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.var_names_make_unique()
    adata.var.index.name = None

    print("Subsetting to oligodendrocytes...")
    mask = adata.obs[CELL_TYPE_COL] == "327 Oligo NN"
    oligo = adata[mask].copy()
    print(f"  {oligo.n_obs} cells, {oligo.n_vars} genes")

    # Preprocess and compute UMAP
    print("Preprocessing...")
    sc.pp.normalize_total(oligo, target_sum=1e4)
    sc.pp.log1p(oligo)
    sc.pp.highly_variable_genes(oligo, n_top_genes=3000, flavor="seurat")
    sc.pp.pca(oligo, n_comps=30, use_highly_variable=True)
    print("Computing KNN + UMAP (this may take a few minutes)...")
    sc.pp.neighbors(oligo, n_neighbors=15, n_pcs=30)
    sc.tl.umap(oligo)
    print("  UMAP done.")

    # UCell scoring on raw counts — reload raw before scoring
    print("Scoring CORT response signature...")
    oligo_raw = adata[mask].copy()
    oligo_raw.var_names = oligo.var_names  # already remapped
    compute_ucell_scores(oligo_raw, SIGNATURE, n_jobs=1)
    oligo.obs["CORT_response_UCell"] = oligo_raw.obs["CORT_response_UCell"].values

    umap = oligo.obsm["X_umap"]
    x, y = umap[:, 0], umap[:, 1]

    # --- Layout: treatment + UCell score + each gene ---
    panel_titles = (
        ["Treatment", "CORT_response (UCell)"]
        + CORT_GENES
    )
    n_panels = len(panel_titles)
    n_cols = 4
    n_rows = int(np.ceil(n_panels / n_cols))

    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=panel_titles,
        horizontal_spacing=0.05,
        vertical_spacing=0.08,
    )

    def rc(i):
        """Convert flat index → (row, col) 1-based."""
        return i // n_cols + 1, i % n_cols + 1

    # Panel 0: treatment
    treatment = oligo.obs[TREATMENT_COL].values
    treatment_colors = {"CORT": "#d62728", "OIL": "#1f77b4"}
    for label, color in treatment_colors.items():
        m = treatment == label
        r, c = rc(0)
        fig.add_trace(go.Scattergl(
            x=x[m], y=y[m], mode="markers",
            marker=dict(size=MARKER_SIZE, opacity=MARKER_OPACITY, color=color),
            name=label, legendgroup="treatment",
            hovertemplate=f"<b>{label}</b><extra></extra>",
        ), row=r, col=c)

    # Panel 1: UCell score
    ucell = oligo.obs["CORT_response_UCell"].values.astype(float)
    r, c = rc(1)
    fig.add_trace(make_colorscale_trace(
        x, y, ucell, "CORT_response (UCell)", colorscale="RdBu_r", zmin=0,
    ), row=r, col=c)

    # Panels 2+: individual genes (log1p expression on normalised data)
    for i, gene in enumerate(CORT_GENES):
        if gene not in oligo.var_names:
            print(f"  WARNING: {gene} not found, skipping")
            continue
        expr = log1p_expr(oligo, gene)
        r, c = rc(i + 2)
        fig.add_trace(make_colorscale_trace(
            x, y, expr, gene, colorscale="Viridis", zmin=0,
        ), row=r, col=c)

    # Style
    fig.update_layout(
        title="Oligodendrocyte UMAP — CORT response signature",
        height=350 * n_rows,
        width=350 * n_cols,
        template="simple_white",
        showlegend=True,
    )
    fig.update_xaxes(showticklabels=False, title_text="")
    fig.update_yaxes(showticklabels=False, title_text="")

    fig.write_html(OUT_HTML, include_plotlyjs="cdn")
    print(f"\nSaved {OUT_HTML}")


if __name__ == "__main__":
    main()
