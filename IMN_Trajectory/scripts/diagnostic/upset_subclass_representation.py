"""
UpSet plots showing subclass representation across conditions for all cells
(no Dcx/Prox1 pre-filtering).

Produces:
  out/subclass_upset_min3_replicates.png
  out/subclass_upset_min2_replicates.png
  out/subclass_upset_min1_replicate.png
  out/subclass_replicate_representation_all_cells.csv

Usage:
    conda run -n anndata_env python scripts/upset_subclass_representation.py
"""

import os
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

# ── config ─────────────────────────────────────────────────────────────────────
DATA          = "data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
OUT_DIR       = "out"
os.makedirs(OUT_DIR, exist_ok=True)

SUBCLASS_COL  = "subclass_name"
SAMPLE_COL    = "sample"
CONDITION_COL = "treatment"
N_REPLICATES  = 3

# ── load ───────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(DATA)
print(f"  Shape: {adata.shape}")

# ── build replicate-representation table ───────────────────────────────────────
print(f"\nBuilding replicate-representation table (condition col: '{CONDITION_COL}')...")

obs = adata.obs[[SUBCLASS_COL, CONDITION_COL, SAMPLE_COL]].copy()
obs[SUBCLASS_COL] = obs[SUBCLASS_COL].astype("category")

conditions     = obs[CONDITION_COL].unique().tolist()
all_subclasses = obs[SUBCLASS_COL].cat.categories.tolist()

rep_counts = (
    obs.groupby([SUBCLASS_COL, CONDITION_COL])[SAMPLE_COL]
    .nunique()
    .unstack(fill_value=0)
    .reindex(index=all_subclasses, fill_value=0)
)


def classify_subclass(row):
    present_in = [c for c in conditions if row[c] >= N_REPLICATES]
    absent_in  = [c for c in conditions if row[c] == 0]
    low_in     = [c for c in conditions if 0 < row[c] < N_REPLICATES]
    if len(present_in) == len(conditions):
        return "keep"
    elif len(present_in) == 1:
        return f"condition_specific ({present_in[0]})"
    elif len(absent_in) > 0:
        return f"absent_in ({', '.join(absent_in)})"
    else:
        return f"low_replicates_in ({', '.join(low_in)})"


rep_counts["status"] = rep_counts.apply(classify_subclass, axis=1)
rep_counts.index.name = SUBCLASS_COL

csv_path = os.path.join(OUT_DIR, "subclass_replicate_representation_all_cells.csv")
rep_counts.to_csv(csv_path)
print(f"  Saved representation table: {csv_path}")

# ── UpSet plots ────────────────────────────────────────────────────────────────
def save_upset(min_reps, out_filename):
    memberships = [
        [c for c in conditions if rep_counts.loc[sc, c] >= min_reps]
        for sc in all_subclasses
    ]
    upset_data = from_memberships(memberships, data=None)
    UpSet(upset_data, subset_size="count", show_counts=True, sort_by="cardinality").plot()
    plt.suptitle(
        f"Subclass representation across conditions (all cells)\n"
        f"Each bar = number of subclasses; intersection = conditions with ≥{min_reps} replicates",
        fontsize=10, y=1.02,
    )
    path = os.path.join(OUT_DIR, out_filename)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")


print("Building UpSet plots...")
save_upset(N_REPLICATES, "subclass_upset_min3_replicates.png")
save_upset(2,            "subclass_upset_min2_replicates.png")
save_upset(1,            "subclass_upset_min1_replicate.png")