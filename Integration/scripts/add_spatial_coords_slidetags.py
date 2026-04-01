"""
add_spatial_coords_slidetags.py
-------------------------------
Add spatial coordinates from per-sample CSV files to the Slide-tags h5ad object.

Join on barcode (cell_bc in CSV, obs index in h5ad — after stripping the '-1' suffix).
Stores x_um and y_um in obs, and a 2-column array in obsm['X_spatial'].
Archives the original h5ad before overwriting.

Usage:
  python add_spatial_coords_slidetags.py
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import shutil

# --- Config ---
DATA_DIR = Path("/scratch/mfafouti/Mommybrain/Integration/data")
H5AD_PATH = DATA_DIR / "Merged_n=207903_all_metadata_slide_tags.h5ad"
ARCHIVE_DIR = DATA_DIR / "archive"

# Map sample names (from obs['sample']) to coordinate CSV files
COORD_FILES = {
    "BC3":  DATA_DIR / "coords_BC3.csv",
    "BC9":  DATA_DIR / "coords_BC9.csv",
    "BC13": DATA_DIR / "coords_BC13.csv",
    "BC14": DATA_DIR / "coords_BC14.csv",
    "BC15": DATA_DIR / "coords_BC15.csv",
    "BC28": DATA_DIR / "coords_BC28.csv",
}

# --- Archive original ---
ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
archive_path = ARCHIVE_DIR / H5AD_PATH.name
if not archive_path.exists():
    print(f"Archiving original to {archive_path}")
    shutil.copy2(H5AD_PATH, archive_path)
else:
    print(f"Archive already exists: {archive_path}")

# --- Load h5ad ---
print(f"Loading {H5AD_PATH}")
adata = sc.read_h5ad(H5AD_PATH)
print(f"Shape: {adata.shape}")

# --- Load and concatenate coordinate CSVs ---
coord_dfs = []
for sample, csv_path in COORD_FILES.items():
    df = pd.read_csv(csv_path)
    df["sample"] = sample
    coord_dfs.append(df)
    print(f"  {sample}: {len(df)} barcodes from {csv_path.name}")

coords = pd.concat(coord_dfs, ignore_index=True)

# Create a barcode column without suffix for matching
# h5ad index: "GGGCATCCACCTCGTG-1" -> strip "-1"
adata.obs["_barcode_bare"] = adata.obs.index.str.replace(r"-\d+$", "", regex=True)

# Join key: barcode + sample (barcodes can repeat across samples)
adata.obs["_join_key"] = adata.obs["_barcode_bare"] + "_" + adata.obs["sample"].astype(str)
coords["_join_key"] = coords["cell_bc"] + "_" + coords["sample"]

# Merge
coords_lookup = coords.set_index("_join_key")[["x_um", "y_um"]]
n_before = len(adata.obs)

adata.obs = adata.obs.join(coords_lookup, on="_join_key", how="left")

# Report match rate
n_matched = adata.obs["x_um"].notna().sum()
print(f"\nMatched {n_matched}/{n_before} cells ({100*n_matched/n_before:.1f}%)")

# Per-sample breakdown
for sample in sorted(COORD_FILES.keys()):
    mask = adata.obs["sample"] == sample
    matched = adata.obs.loc[mask, "x_um"].notna().sum()
    total = mask.sum()
    print(f"  {sample}: {matched}/{total} ({100*matched/total:.1f}%)")

# Check for cells with (0, 0) coordinates — these are unlocalized in the CSVs
zero_mask = (adata.obs["x_um"] == 0.0) & (adata.obs["y_um"] == 0.0)
n_zero = zero_mask.sum()
if n_zero > 0:
    print(f"\nNote: {n_zero} cells have (0, 0) coordinates (unlocalized) — setting to NaN")
    adata.obs.loc[zero_mask, ["x_um", "y_um"]] = np.nan

# Store in obsm as well
spatial = adata.obs[["x_um", "y_um"]].values.astype(np.float32)
adata.obsm["X_spatial"] = spatial

# Clean up temp columns
adata.obs.drop(columns=["_barcode_bare", "_join_key"], inplace=True)

# --- Save ---
print(f"\nSaving updated object to {H5AD_PATH}")
adata.write_h5ad(H5AD_PATH)
print("Done.")
