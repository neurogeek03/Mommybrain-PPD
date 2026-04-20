"""
Export spatial viewer data as compact binary files + viewer.html.
Usage: python export_binary.py --config test.conf

Output layout under out_dir/:
  manifest.json           sample metadata, gene list, cell type palette
  {sample}.bin            per-sample: x (float32), y (float32), celltype_id (uint16)
                          layout: [x×n][y×n][ct_id×n] = 10×n bytes
  genes/{gene}.bin        per-gene: log1p-normalized expression (float32)
                          layout: all samples concatenated in config order
  viewer.html             standalone viewer (fetches the above files)
"""

import argparse
import ast
import configparser
import json
import os
import shutil

import numpy as np
import scipy.sparse
from plotly.colors import qualitative

from utils.data_utils import (
    load_adata,
    normalize,
    subset_samples,
    get_celltypes,
    load_gene_list,
    get_all_gene_symbols,
    load_celltype_colors,
)


# ── Config ─────────────────────────────────────────────────────────────────────

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
        "out_dir":      cfg.get("out_dir", "viewer_data"),
    }


parser = argparse.ArgumentParser(description="Export spatial viewer binary data")
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg  = parse_config(args.config)

SAMPLES  = cfg["samples"]
out_dir  = cfg["out_dir"]
gene_dir = os.path.join(out_dir, "genes")
os.makedirs(gene_dir, exist_ok=True)


# ── Load & preprocess ──────────────────────────────────────────────────────────

print(f"Loading {cfg['h5ad_path']} ...")
adata = load_adata(cfg["h5ad_path"])

print(f"Subsetting to samples: {SAMPLES} ...")
adata = subset_samples(adata, SAMPLES, cfg["sample_col"])

print("Normalizing ...")
adata = normalize(adata)

print("Building cell type palette ...")
celltypes = get_celltypes(adata, cfg["celltype_col"])
if cfg["color_csv"]:
    color_map = load_celltype_colors(cfg["color_csv"])
    n_matched = sum(1 for ct in celltypes if ct in color_map)
    print(f"  Loaded colors: {n_matched}/{len(celltypes)} matched")
else:
    color_pool = qualitative.Alphabet + qualitative.Dark24
    color_map  = {ct: color_pool[i % len(color_pool)] for i, ct in enumerate(celltypes)}

celltype_to_id  = {ct: i for i, ct in enumerate(celltypes)}
celltype_palette = [{"name": ct, "color": color_map.get(ct, "#808080")} for ct in celltypes]
assert len(celltypes) < 2**16, "Too many cell types for uint16 encoding"

if cfg["gene_list"]:
    print(f"Loading gene list from {cfg['gene_list']} ...")
    genes = load_gene_list(cfg["gene_list"], adata)
    print(f"  {len(genes)} genes loaded")
else:
    genes = get_all_gene_symbols(adata)


# ── Fast lookup structures ─────────────────────────────────────────────────────

print("Converting X to CSC for fast gene extraction ...")
if scipy.sparse.issparse(adata.X):
    adata.X = scipy.sparse.csc_matrix(adata.X)

symbol_to_idx  = {sym: i for i, sym in enumerate(adata.var["gene_symbol"].values)}
obs_sample     = adata.obs[cfg["sample_col"]].values
sample_indices = {s: np.where(obs_sample == s)[0] for s in SAMPLES}
coords         = adata.obsm[cfg["spatial_key"]]
celltype_arr   = adata.obs[cfg["celltype_col"]].astype(str).values


# ── Per-sample binary files ────────────────────────────────────────────────────
# Layout: [x × n float32][y × n float32][celltype_id × n uint16] = 10n bytes

print("Writing per-sample binary files ...")
sample_meta = []
for sample in SAMPLES:
    idx = sample_indices[sample]
    n   = len(idx)

    x      = coords[idx, 0].astype(np.float32)
    y      = coords[idx, 1].astype(np.float32)
    ct_ids = np.array([celltype_to_id[celltype_arr[i]] for i in idx], dtype=np.uint16)

    buf  = bytearray(10 * n)
    view = memoryview(buf)
    np.frombuffer(view[0      : 4*n],  dtype=np.float32)[:] = x
    np.frombuffer(view[4*n    : 8*n],  dtype=np.float32)[:] = y
    np.frombuffer(view[8*n    : 10*n], dtype=np.uint16)[:] = ct_ids

    out_path = os.path.join(out_dir, f"{sample}.bin")
    with open(out_path, "wb") as f:
        f.write(buf)

    sample_meta.append({
        "name":    sample,
        "n_cells": n,
        "x_range": [float(x.min()), float(x.max())],
        "y_range": [float(y.min()), float(y.max())],
        "file":    f"{sample}.bin",
    })
    print(f"  {sample}: {n:,} cells → {sample}.bin  ({10*n/1e6:.1f} MB)")


# ── Per-gene binary files ──────────────────────────────────────────────────────
# Layout: float32 expression values for all samples concatenated in SAMPLES order

print(f"Writing {len(genes)} gene binary files ...")
for j, gene in enumerate(genes):
    col_idx = symbol_to_idx[gene]
    col = adata.X[:, col_idx]
    if scipy.sparse.issparse(col):
        col = np.asarray(col.todense()).ravel()
    else:
        col = np.asarray(col).ravel()

    expr = np.concatenate([col[sample_indices[s]] for s in SAMPLES]).astype(np.float32)
    with open(os.path.join(gene_dir, f"{gene}.bin"), "wb") as f:
        f.write(expr.tobytes())

    if (j + 1) % 100 == 0 or j + 1 == len(genes):
        print(f"  {j+1}/{len(genes)}")


# ── manifest.json ──────────────────────────────────────────────────────────────

manifest = {
    "version":          1,
    "point_size":       cfg["point_size"],
    "samples":          sample_meta,
    "celltype_palette": celltype_palette,
    "genes":            genes,
}

manifest_path = os.path.join(out_dir, "manifest.json")
with open(manifest_path, "w") as f:
    json.dump(manifest, f, indent=2)
print(f"manifest.json written ({len(genes)} genes, {len(celltypes)} cell types)")


# ── Copy viewer.html ───────────────────────────────────────────────────────────

viewer_src = os.path.join(os.path.dirname(__file__), "viewer.html")
viewer_dst = os.path.join(out_dir, "viewer.html")
shutil.copy(viewer_src, viewer_dst)
print(f"viewer.html copied to {viewer_dst}")

total_mb = sum(
    os.path.getsize(os.path.join(dp, fn))
    for dp, _, fns in os.walk(out_dir)
    for fn in fns
) / 1e6
print(f"\nDone. Total output size: {total_mb:.0f} MB in {out_dir}/")