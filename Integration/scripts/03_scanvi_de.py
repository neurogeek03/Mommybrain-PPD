"""
03_scanvi_de.py
---------------
scANVI posterior-based differential expression.
Primary comparison: CORT vs OIL at PD23, per cell type.

Loads the saved scANVI model and integrated.h5ad.
For each cell type with sufficient cells in both groups, runs
model.differential_expression() using boolean index masks.
Batch correction is implicit — the model was trained with sample_id as batch_key.

Usage:
  python 03_scanvi_de.py [config_path]
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import yaml
from pathlib import Path

# =============================================================================
# CONFIG
# =============================================================================

parser = argparse.ArgumentParser()
parser.add_argument(
    "config",
    nargs="?",
    default=Path(__file__).parent.parent / "config" / "03_scanvi_de_configB.conf",
)
args = parser.parse_args()

with open(args.config) as f:
    cfg = yaml.safe_load(f)

INTEGRATED_PATH    = cfg["integrated_path"]
MODEL_PATH         = cfg["scanvi_model_path"]
OUT_DIR            = Path(cfg["out_dir"])
CELL_TYPE_KEY      = cfg["cell_type_key"]
TREATMENT_KEY      = cfg["treatment_key"]
DAY_KEY            = cfg["day_key"]
DAY_FILTER         = cfg["day_filter"]
GROUP1             = cfg["group1"]
GROUP2             = cfg["group2"]
MIN_CELLS          = cfg["min_cells_per_group"]
REQUIRE_BOTH_METHODS  = cfg.get("require_both_methods", False)
MAX_CELLS          = cfg.get("max_cells_per_group")     # may be None
N_SAMPLES          = cfg.get("n_samples", 1000)
BATCH_SIZE         = cfg.get("batch_size", 256)
LFC_THRESH         = cfg.get("lfc_threshold", 0.25)
PROBA_THRESH       = cfg.get("proba_de_threshold", 0.8)

OUT_DIR.mkdir(parents=True, exist_ok=True)
per_ct_dir = OUT_DIR / "per_cell_type"
per_ct_dir.mkdir(exist_ok=True)

# =============================================================================
# LOAD DATA + MODEL
# =============================================================================

print(f"\nLoading integrated object: {INTEGRATED_PATH}")
adata = sc.read_h5ad(INTEGRATED_PATH)
print(f"  Shape: {adata.shape}")

print(f"\nLoading scANVI model: {MODEL_PATH}")
model = scvi.model.SCANVI.load(MODEL_PATH, adata=adata)
print("  Model loaded.")

# =============================================================================
# DEFINE CELL GROUPS
# =============================================================================

# Restrict to the target day
day_mask = adata.obs[DAY_KEY] == DAY_FILTER
print(f"\nCells at {DAY_FILTER}: {day_mask.sum():,} / {len(adata):,}")

g1_mask = day_mask & (adata.obs[TREATMENT_KEY] == GROUP1)
g2_mask = day_mask & (adata.obs[TREATMENT_KEY] == GROUP2)
print(f"  {GROUP1}: {g1_mask.sum():,}   {GROUP2}: {g2_mask.sum():,}")

cell_types = sorted(adata.obs[CELL_TYPE_KEY].unique())
print(f"\nCell types in object: {len(cell_types)}")

# =============================================================================
# PER-CELL-TYPE DE
# =============================================================================

stats_rows = []
skipped = []

for ct in cell_types:
    ct_mask = adata.obs[CELL_TYPE_KEY] == ct

    idx1 = np.where(g1_mask & ct_mask)[0]
    idx2 = np.where(g2_mask & ct_mask)[0]

    n1, n2 = len(idx1), len(idx2)

    if n1 < MIN_CELLS or n2 < MIN_CELLS:
        reason = f"n_{GROUP1}={n1}, n_{GROUP2}={n2} (min={MIN_CELLS})"
        print(f"  SKIP  {ct:50s}  {reason}")
        skipped.append({"cell_type": ct, "n_group1": n1, "n_group2": n2, "reason": reason})
        continue

    # Both-methods filter: cell type must appear in slide_tags AND slide_seq
    if REQUIRE_BOTH_METHODS:
        methods = set(adata.obs.loc[(g1_mask | g2_mask) & ct_mask, "method"].unique())
        if not {"slide_tags", "slide_seq"}.issubset(methods):
            reason = f"only present in: {methods}"
            print(f"  SKIP  {ct:50s}  {reason}")
            skipped.append({"cell_type": ct, "n_group1": n1, "n_group2": n2, "reason": reason})
            continue

    # Optional: subsample to cap memory / runtime
    rng = np.random.default_rng(42)
    if MAX_CELLS and n1 > MAX_CELLS:
        idx1 = rng.choice(idx1, MAX_CELLS, replace=False)
    if MAX_CELLS and n2 > MAX_CELLS:
        idx2 = rng.choice(idx2, MAX_CELLS, replace=False)

    print(f"  DE    {ct:50s}  n1={len(idx1):5d}  n2={len(idx2):5d}", end="  ", flush=True)

    de_df = model.differential_expression(
        idx1=idx1,
        idx2=idx2,
        batch_size=BATCH_SIZE,
    )

    # Compute LFC from scale1/scale2 (normalized expression per group)
    # proba_m1 = P(higher in group1/CORT), proba_m2 = P(higher in group2/OIL)
    eps = 1e-10
    de_df["lfc"] = np.log2(de_df["scale1"] + eps) - np.log2(de_df["scale2"] + eps)

    # proba_de: probability of being DE in either direction
    de_df["proba_de"] = de_df[["proba_m1", "proba_m2"]].max(axis=1)

    # Save full results for this cell type
    safe_name = ct.replace("/", "-").replace(" ", "_")
    out_path = per_ct_dir / f"de_{safe_name}.csv"
    de_df.to_csv(out_path)

    # Summarize
    n_de = int((de_df["proba_de"] > PROBA_THRESH).sum())
    n_de_lfc = int(
        ((de_df["proba_de"] > PROBA_THRESH) & (de_df["lfc"].abs() > LFC_THRESH)).sum()
    )
    median_lfc = de_df["lfc"].median()
    print(f"n_de(p>{PROBA_THRESH})={n_de}  n_de+lfc={n_de_lfc}  median_lfc={median_lfc:.3f}")

    stats_rows.append({
        "cell_type":      ct,
        "n_group1":       n1,
        "n_group2":       n2,
        "n_used_g1":      len(idx1),
        "n_used_g2":      len(idx2),
        "n_genes_tested": len(de_df),
        f"n_de_p{PROBA_THRESH}":              n_de,
        f"n_de_p{PROBA_THRESH}_lfc{LFC_THRESH}": n_de_lfc,
        "median_lfc":     median_lfc,
    })

# =============================================================================
# SUMMARY OUTPUTS
# =============================================================================

stats_df = pd.DataFrame(stats_rows)
stats_path = OUT_DIR / "de_stats.csv"
stats_df.to_csv(stats_path, index=False)
print(f"\nDE stats saved: {stats_path}")
print(stats_df.to_string(index=False))

if skipped:
    skip_df = pd.DataFrame(skipped)
    skip_path = OUT_DIR / "de_skipped.csv"
    skip_df.to_csv(skip_path, index=False)
    print(f"\nSkipped cell types ({len(skipped)}): {skip_path}")

print(f"\n[03_scanvi_de.py] Done. Outputs in: {OUT_DIR}")
