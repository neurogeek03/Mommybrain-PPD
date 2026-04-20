"""
Spatial Gene Expression Interactive Viewer — Dash app
Usage: python main.py --config test.conf
Then open http://localhost:<port> in your browser.
"""

import argparse
import ast
import configparser

import numpy as np
from plotly.colors import qualitative

import scipy.sparse

import dash
from dash import dcc, html, Input, Output, Patch

from utils.data_utils import (
    load_adata,
    normalize,
    subset_samples,
    get_spatial_per_sample,
    get_celltype_per_sample,
    get_celltypes,
    get_all_gene_symbols,
    load_gene_list,
    load_celltype_colors,
)
from utils.plot_utils import build_figure


# ── Config ────────────────────────────────────────────────────────────────────

def parse_config(path):
    config = configparser.ConfigParser()
    config.read(path)
    cfg = config["settings"]
    return {
        "h5ad_path":    cfg["h5ad_path"],
        "samples":      ast.literal_eval(cfg["samples"]),
        "sample_col":   cfg.get("sample_col", "sample"),
        "celltype_col": cfg.get("celltype_col", "RCTD_first_type_rat"),
        "spatial_key":  cfg.get("spatial_key", "X_spatial"),
        "point_size":   int(cfg.get("point_size", "2")),
        "port":         int(cfg.get("port", "8050")),
        "debug":        cfg.get("debug", "false").lower() == "true",
        "color_csv":    cfg.get("color_csv", ""),
        "gene_list":    cfg.get("gene_list", ""),
    }


parser = argparse.ArgumentParser(description="Spatial gene expression viewer")
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg = parse_config(args.config)


# ── Load & preprocess data (once at startup) ──────────────────────────────────

print(f"Loading {cfg['h5ad_path']} ...")
adata = load_adata(cfg["h5ad_path"])

print(f"Subsetting to samples: {cfg['samples']} ...")
adata = subset_samples(adata, cfg["samples"], cfg["sample_col"])

print("Normalizing (normalize_total + log1p) ...")
adata = normalize(adata)

print("Pre-computing spatial coordinates and cell type labels ...")
spatial_data   = get_spatial_per_sample(adata, cfg["samples"], cfg["sample_col"], cfg["spatial_key"])
celltype_data  = get_celltype_per_sample(adata, cfg["samples"], cfg["sample_col"], cfg["celltype_col"])
celltypes      = get_celltypes(adata, cfg["celltype_col"])
if cfg["gene_list"]:
    print(f"Loading gene list from {cfg['gene_list']} ...")
    all_genes = load_gene_list(cfg["gene_list"], adata)
    print(f"  {len(all_genes)} genes loaded from list")
else:
    all_genes = get_all_gene_symbols(adata)

# Categorical color map — load from CSV if provided, else cycle qualitative palettes
if cfg["color_csv"]:
    celltype_color_map = load_celltype_colors(cfg["color_csv"])
    n_matched = sum(1 for ct in celltypes if ct in celltype_color_map)
    print(f"  Loaded colors from CSV: {n_matched}/{len(celltypes)} cell types matched")
else:
    color_pool = qualitative.Alphabet + qualitative.Dark24
    celltype_color_map = {ct: color_pool[i % len(color_pool)] for i, ct in enumerate(celltypes)}

SAMPLES    = cfg["samples"]
POINT_SIZE = cfg["point_size"]

# ── Performance: pre-compute structures used in every callback ─────────────────

# 1. CSC format — makes column (gene) extraction fast instead of scanning all rows
print("Converting X to CSC format for fast gene lookup ...")
if scipy.sparse.issparse(adata.X):
    adata.X = scipy.sparse.csc_matrix(adata.X)

# 2. Pre-computed row indices per sample — avoids rebuilding masks on every callback
obs_sample = adata.obs[cfg["sample_col"]].values
SAMPLE_INDICES = {
    sample: np.where(obs_sample == sample)[0]
    for sample in SAMPLES
}

# 3. Pre-built gene symbol → column index map
SYMBOL_TO_IDX = {
    sym: i for i, sym in enumerate(adata.var["gene_symbol"].values)
}

print(f"Ready. {len(all_genes)} genes available, {len(celltypes)} cell types.")


# ── Helper: fast single-gene expression extraction ────────────────────────────

def get_gene_expr_fast(gene):
    """Extract expression for one gene across all samples using pre-built structures."""
    col_idx = SYMBOL_TO_IDX[gene]
    col = adata.X[:, col_idx]
    if scipy.sparse.issparse(col):
        col = np.asarray(col.todense()).ravel()
    else:
        col = np.asarray(col).ravel()
    return {sample: col[indices].tolist() for sample, indices in SAMPLE_INDICES.items()}


# ── Helper: compute per-sample color arrays ───────────────────────────────────

def compute_colors(gene, colormode, filter_ct):
    """
    Returns (color_data, colorscale, cmin, cmax).
    color_data : {sample: list} — floats or hex strings; None = transparent.
    """
    if colormode == "expression":
        expr_by_sample = get_gene_expr_fast(gene)
        all_vals = [v for vals in expr_by_sample.values() for v in vals]
        cmin = float(np.min(all_vals)) if all_vals else 0.0
        cmax = float(np.max(all_vals)) if all_vals else 1.0

        color_data = {}
        for sample in SAMPLES:
            vals = expr_by_sample[sample]
            if filter_ct == "All":
                color_data[sample] = vals
            else:
                cts = celltype_data[sample]
                color_data[sample] = [v if ct == filter_ct else None
                                       for v, ct in zip(vals, cts)]
        return color_data, "Viridis", cmin, cmax

    else:  # celltype
        color_data = {}
        for sample in SAMPLES:
            cts = celltype_data[sample]
            if filter_ct == "All":
                color_data[sample] = [celltype_color_map.get(ct, "#888888") for ct in cts]
            else:
                color_data[sample] = [
                    celltype_color_map.get(ct, "#888888") if ct == filter_ct else None
                    for ct in cts
                ]
        return color_data, None, None, None


# ── Initial figure ────────────────────────────────────────────────────────────

initial_gene = all_genes[0]
init_colors, init_cs, init_cmin, init_cmax = compute_colors(initial_gene, "expression", "All")
initial_fig = build_figure(
    SAMPLES, spatial_data, init_colors, celltype_data,
    colorscale=init_cs, cmin=init_cmin, cmax=init_cmax,
    point_size=POINT_SIZE,
)


# ── Dash layout ───────────────────────────────────────────────────────────────

app = dash.Dash(__name__, title="Spatial Viewer")

app.layout = html.Div(
    style={"background": "#1a1a1a", "minHeight": "100vh", "fontFamily": "'Segoe UI', sans-serif"},
    children=[
        # Controls bar
        html.Div(
            style={
                "display": "flex", "alignItems": "center", "gap": "24px",
                "padding": "12px 20px", "background": "#2a2a2a",
                "borderBottom": "1px solid #444", "flexWrap": "wrap",
            },
            children=[
                html.Div([
                    html.Label("Gene", style={"color": "#aaa", "fontSize": "12px",
                                              "textTransform": "uppercase", "letterSpacing": "0.05em"}),
                    dcc.Dropdown(
                        id="gene-select",
                        options=[{"label": g, "value": g} for g in all_genes],
                        value=initial_gene,
                        clearable=False,
                        searchable=True,
                        style={"minWidth": "220px", "background": "#333", "color": "#000"},
                    ),
                ], style={"display": "flex", "flexDirection": "column", "gap": "4px"}),

                html.Div([
                    html.Label("Color by", style={"color": "#aaa", "fontSize": "12px",
                                                  "textTransform": "uppercase", "letterSpacing": "0.05em"}),
                    dcc.RadioItems(
                        id="colormode",
                        options=[
                            {"label": " Expression", "value": "expression"},
                            {"label": " Cell Type",  "value": "celltype"},
                        ],
                        value="expression",
                        inline=True,
                        style={"color": "#e0e0e0", "fontSize": "13px"},
                        inputStyle={"marginRight": "4px", "accentColor": "#7ec8e3"},
                        labelStyle={"marginRight": "14px"},
                    ),
                ], style={"display": "flex", "flexDirection": "column", "gap": "4px"}),

                html.Div([
                    html.Label("Filter cell type", style={"color": "#aaa", "fontSize": "12px",
                                                          "textTransform": "uppercase", "letterSpacing": "0.05em"}),
                    dcc.Dropdown(
                        id="celltype-select",
                        options=[{"label": "All cell types", "value": "All"}]
                                + [{"label": ct, "value": ct} for ct in celltypes],
                        value="All",
                        clearable=False,
                        searchable=True,
                        style={"minWidth": "260px", "background": "#333", "color": "#000"},
                    ),
                ], style={"display": "flex", "flexDirection": "column", "gap": "4px"}),
            ],
        ),

        # Plot
        dcc.Graph(
            id="spatial-plot",
            figure=initial_fig,
            style={"width": "100%"},
            config={"displayModeBar": True, "scrollZoom": True},
        ),
    ],
)


# ── Callback ──────────────────────────────────────────────────────────────────

@app.callback(
    Output("spatial-plot", "figure"),
    Input("gene-select",    "value"),
    Input("colormode",      "value"),
    Input("celltype-select","value"),
)
def update_plot(gene, colormode, filter_ct):
    colors, colorscale, cmin, cmax = compute_colors(gene, colormode, filter_ct)
    # Patch sends only the changed marker properties — not the full figure with x/y coords
    patched = Patch()
    for i, sample in enumerate(SAMPLES):
        is_last = (i == len(SAMPLES) - 1)
        patched["data"][i]["marker"]["color"]      = colors[sample]
        patched["data"][i]["marker"]["colorscale"] = colorscale
        patched["data"][i]["marker"]["cmin"]       = cmin
        patched["data"][i]["marker"]["cmax"]       = cmax
        patched["data"][i]["marker"]["showscale"]  = is_last and colormode == "expression"
    return patched


# ── Run ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"Starting Dash app on http://localhost:{cfg['port']}")
    app.run(host="0.0.0.0", port=cfg["port"], debug=cfg["debug"])
