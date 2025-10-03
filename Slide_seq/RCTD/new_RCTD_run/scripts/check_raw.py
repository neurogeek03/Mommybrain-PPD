#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scanpy as sc
import numpy as np
import scipy.sparse as sp
import argparse

parser = argparse.ArgumentParser(description="Check if AnnData contains raw counts")
parser.add_argument("--input", required=True, help="Input .h5ad file")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

# choose which matrix to inspect
if adata.raw is not None:
    X = adata.raw.X
    print("Using adata.raw.X as counts")
elif "counts" in adata.layers:
    X = adata.layers["counts"]
    print("Using adata.layers['counts'] as counts")
else:
    X = adata.X
    print("WARNING: using adata.X (may be normalized)")

# convert sparse to dense if needed
if sp.issparse(X):
    X = X.toarray()

# check min, max, integer status
print(f"Shape: {X.shape}")
print(f"Min value: {X.min()}")
print(f"Max value: {X.max()}")
print(f"Are all values integers? {np.all(np.mod(X,1)==0)}")
