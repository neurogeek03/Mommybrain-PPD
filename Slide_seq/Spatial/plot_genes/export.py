"""
Export standalone HTML spatial gene expression viewer.
Usage: python export.py --config test.conf
"""

import argparse
import ast
import configparser

import numpy as np
import scipy.sparse
from plotly.colors import qualitative

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
from utils.html_utils import build_html


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
        "color_csv":    cfg.get("color_csv", ""),
        "gene_list":    cfg.get("gene_list", ""),
        "output_html":  cfg["output_html"],
    }


parser = argparse.ArgumentParser(description="Export spatial viewer to standalone HTML")
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg = parse_config(args.config)

SAMPLES    = cfg["samples"]
POINT_SIZE = cfg["point_size"]


# ── Load & preprocess ─────────────────────────────────────────────────────────

print(f"Loading {cfg['h5ad_path']} ...")
adata = load_adata(cfg["h5ad_path"])

print(f"Subsetting to samples: {SAMPLES} ...")
adata = subset_samples(adata, SAMPLES, cfg["sample_col"])

print("Normalizing ...")
adata = normalize(adata)

print("Pre-computing spatial coordinates and cell type labels ...")
spatial_data  = get_spatial_per_sample(adata, SAMPLES, cfg["sample_col"], cfg["spatial_key"])
celltype_data = get_celltype_per_sample(adata, SAMPLES, cfg["sample_col"], cfg["celltype_col"])
celltypes     = get_celltypes(adata, cfg["celltype_col"])

if cfg["gene_list"]:
    print(f"Loading gene list from {cfg['gene_list']} ...")
    all_genes = load_gene_list(cfg["gene_list"], adata)
    print(f"  {len(all_genes)} genes loaded")
else:
    all_genes = get_all_gene_symbols(adata)

if cfg["color_csv"]:
    celltype_color_map = load_celltype_colors(cfg["color_csv"])
    n_matched = sum(1 for ct in celltypes if ct in celltype_color_map)
    print(f"Loaded colors from CSV: {n_matched}/{len(celltypes)} cell types matched")
else:
    color_pool = qualitative.Alphabet + qualitative.Dark24
    celltype_color_map = {ct: color_pool[i % len(color_pool)] for i, ct in enumerate(celltypes)}


# ── Pre-compute expression for all genes ─────────────────────────────────────

print("Converting X to CSC for fast gene lookup ...")
if scipy.sparse.issparse(adata.X):
    adata.X = scipy.sparse.csc_matrix(adata.X)

symbol_to_idx = {sym: i for i, sym in enumerate(adata.var["gene_symbol"].values)}
obs_sample    = adata.obs[cfg["sample_col"]].values
sample_indices = {s: np.where(obs_sample == s)[0] for s in SAMPLES}

print(f"Extracting expression for {len(all_genes)} genes ...")
expr_data = {}
for gene in all_genes:
    col_idx = symbol_to_idx[gene]
    col = adata.X[:, col_idx]
    if scipy.sparse.issparse(col):
        col = np.asarray(col.todense()).ravel()
    else:
        col = np.asarray(col).ravel()
    expr_data[gene] = {s: col[sample_indices[s]].tolist() for s in SAMPLES}

all_vals   = [v for gd in expr_data.values() for sd in gd.values() for v in sd]
global_min = float(np.min(all_vals))
global_max = float(np.max(all_vals))
print(f"Expression range: [{global_min:.3f}, {global_max:.3f}]")


# ── Build initial figure (first gene, expression mode) ────────────────────────

initial_gene = all_genes[0]
init_colors  = {s: expr_data[initial_gene][s] for s in SAMPLES}

initial_fig = build_figure(
    SAMPLES, spatial_data, init_colors, celltype_data,
    colorscale="Viridis", cmin=global_min, cmax=global_max,
    point_size=POINT_SIZE,
)


# ── Write HTML ────────────────────────────────────────────────────────────────

print(f"Writing HTML to {cfg['output_html']} ...")
build_html(
    fig=initial_fig,
    expr_data=expr_data,
    celltype_data=celltype_data,
    genes=all_genes,
    celltypes=celltypes,
    sample_names=SAMPLES,
    output_path=cfg["output_html"],
    point_size=POINT_SIZE,
    global_min=global_min,
    global_max=global_max,
    celltype_color_map=celltype_color_map,
)
