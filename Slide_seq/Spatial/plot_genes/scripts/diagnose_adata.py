"""
Diagnostic script — inspect h5ad structure before building the viewer.
Run with: python diagnose_adata.py
"""
import scanpy as sc
import numpy as np

path = "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/data/All_RCTD_types_singlet_score_0_slide_seq_15.h5ad"

print("Loading...")
adata = sc.read_h5ad(path)

print("\n=== AnnData summary ===")
print(adata)

print("\n=== obs columns ===")
print(adata.obs.columns.tolist())

print("\n=== obsm keys ===")
print(list(adata.obsm.keys()))

print("\n=== layers ===")
print(list(adata.layers.keys()))

print("\n=== var columns ===")
print(adata.var.columns.tolist())
print("var_names sample:", adata.var_names[:5].tolist())

print("\n=== Sample column check ===")
if "sample" in adata.obs.columns:
    samples = adata.obs["sample"].unique()
    print(f"Found {len(samples)} samples:")
    print(sorted(samples))
else:
    print("No 'sample' column — check obs columns above")

print("\n=== Cell type column check ===")
ct_col = "RCTD_first_type_rat"
if ct_col in adata.obs.columns:
    cts = adata.obs[ct_col].unique()
    print(f"Found {len(cts)} cell types in '{ct_col}':")
    print(sorted(cts))
else:
    print(f"No '{ct_col}' column — check obs columns above")

print("\n=== Spatial key check ===")
for key in adata.obsm.keys():
    arr = adata.obsm[key]
    print(f"  {key}: shape {arr.shape}, dtype {arr.dtype}")

print("\n=== X matrix ===")
print(f"  dtype: {adata.X.dtype}")
print(f"  shape: {adata.X.shape}")
print(f"  min: {adata.X.min():.4f}, max: {adata.X.max():.4f}")
# Check if likely raw (integers) or normalized (floats)
import scipy.sparse as sp
sample_vals = adata.X[:100, :100]
if sp.issparse(sample_vals):
    sample_vals = sample_vals.toarray()
print(f"  Sample values (first 10 non-zero): {sample_vals[sample_vals > 0].flat[:10]}")

print("\n=== var 'name' column check ===")
if "name" in adata.var.columns:
    name_col = adata.var["name"].dropna().astype(str)
    print(f"  First 10 values: {name_col[:10].tolist()}")
    # Check for dot-separated format gene.ENSEMBL
    dot_pattern = name_col.str.match(r"^.+\.[A-Z0-9]+$")
    print(f"  Entries matching 'name.ID' pattern: {dot_pattern.sum()} / {len(name_col)}")
    if dot_pattern.any():
        examples = name_col[dot_pattern][:5].tolist()
        print(f"  Examples: {examples}")
        # Preview what the gene name part would be
        extracted = name_col[dot_pattern].str.rsplit(".", n=1).str[0]
        print(f"  Extracted name part (before last dot): {extracted[:5].tolist()}")
else:
    print("  No 'name' column found")

print("\n=== Gene symbol check ===")
if "gene_symbol" in adata.var.columns:
    syms = adata.var["gene_symbol"].dropna().astype(str)
    proper = syms[~syms.str.match(r"^ENS[A-Z]+\d+$")]
    ensembl_only = syms[syms.str.match(r"^ENS[A-Z]+\d+$")]
    print(f"  Total genes:               {len(adata.var)}")
    print(f"  Proper gene symbols:       {proper.nunique()}")
    print(f"  Ensembl ID placeholders:   {ensembl_only.nunique()}")
    print(f"  Examples (proper):         {proper.unique()[:10].tolist()}")
    print(f"  Examples (Ensembl-only):   {ensembl_only.unique()[:5].tolist()}")
else:
    print("  No 'gene_symbol' column found in adata.var")