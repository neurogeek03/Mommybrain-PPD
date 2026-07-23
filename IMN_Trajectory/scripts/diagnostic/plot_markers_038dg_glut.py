"""
Bar chart of Dcx, Prox1, Calb1, Calb2 expression in '038 DG Glut' cells.

Shows, per sample, how many cells express each gene individually
and all three together (once with Calb1, once with Calb2).

Usage:
    conda run -n anndata_env python scripts/diagnostic/plot_markers_038dg_glut.py \
        data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad
"""

import argparse
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── config ─────────────────────────────────────────────────────────────────────
SUBCLASS_COL  = "subclass_name"
SAMPLE_COL    = "sample"
TARGET_CLASS  = "037 DG Glut"   # partial match used if exact not found

GENES         = ["Dcx", "Prox1", "Calb1", "Calb2"]

CORT_SAMPLES  = {"BC28", "BC3", "BC9"}
OIL_SAMPLES   = {"BC15", "BC14", "BC13"}
SAMPLE_ORDER  = ["BC13", "BC14", "BC15", "BC3", "BC9", "BC28"]

COLORS = {
    "Dcx":               "#2980B9",
    "Prox1":             "#8E44AD",
    "Calb1":             "#E67E22",
    "Calb2":             "#27AE60",
    "Dcx+Prox1+Calb1":  "#C0392B",
    "Dcx+Prox1+Calb2":  "#922B21",
}

OUT_DIR  = "out"
OUT_FILE = "markers_037dg_glut_dcx_prox1_calb1_calb2.png"

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("h5ad", nargs="?",
    default="data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad")
args = parser.parse_args()

# ── load ───────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(args.h5ad)
print(f"  Shape: {adata.shape}")

# ── subset to target subclass (exact first, then partial) ─────────────────────
labels = adata.obs[SUBCLASS_COL].astype(str)
exact  = labels == TARGET_CLASS
if exact.sum() == 0:
    partial = labels.str.contains(TARGET_CLASS, case=False, na=False)
    if partial.sum() == 0:
        available = sorted(labels.unique())
        raise ValueError(
            f"No cells match '{TARGET_CLASS}'. Available subclasses:\n  "
            + "\n  ".join(available)
        )
    matched = labels[partial].unique()
    if len(matched) > 1:
        raise ValueError(f"Ambiguous partial match for '{TARGET_CLASS}': {matched.tolist()}")
    found = matched[0]
    print(f"  Partial match: '{TARGET_CLASS}' → '{found}'")
    mask_class = partial
else:
    found = TARGET_CLASS
    mask_class = exact

sub = adata[mask_class].copy()
print(f"  '{found}' cells: {sub.n_obs}")

# ── locate gene_symbol column ──────────────────────────────────────────────────
sym_cols = [c for c in sub.var.columns if "gene_symbol" in c.lower()]
if not sym_cols:
    raise ValueError(f"No gene_symbol column in adata.var. Available: {sub.var.columns.tolist()}")
sym_col = sym_cols[0]
symbols = sub.var[sym_col].astype(str)


def resolve_gene(name):
    exact = symbols.str.lower() == name.lower()
    hits  = sub.var_names[exact].tolist()
    if not hits:
        partial = symbols.str.lower().str.contains(name.lower(), na=False)
        hits    = sub.var_names[partial].tolist()
    if not hits:
        raise ValueError(f"Gene '{name}' not found in adata.var['{sym_col}'].")
    if len(hits) > 1:
        print(f"  WARNING: '{name}' matched {len(hits)} entries: {hits} — using first.")
    label = symbols[hits[0]]
    print(f"  Resolved '{name}' → var_name: {hits[0]}  symbol: {label}")
    return hits[0]


def get_expr(var_name):
    idx = list(sub.var_names).index(var_name)
    col = sub.X[:, idx]
    return np.asarray(col.todense()).flatten() if sp.issparse(col) else np.asarray(col).flatten()


# ── resolve genes and build per-cell boolean masks ────────────────────────────
print("\nResolving genes...")
gene_masks = {}
for g in GENES:
    vname = resolve_gene(g)
    expr  = get_expr(vname)
    gene_masks[g] = expr > 0
    print(f"  {g}+: {gene_masks[g].sum()} / {sub.n_obs} cells")

all_calb1 = gene_masks["Dcx"] & gene_masks["Prox1"] & gene_masks["Calb1"]
all_calb2 = gene_masks["Dcx"] & gene_masks["Prox1"] & gene_masks["Calb2"]
print(f"  Dcx+Prox1+Calb1+: {all_calb1.sum()} / {sub.n_obs} cells")
print(f"  Dcx+Prox1+Calb2+: {all_calb2.sum()} / {sub.n_obs} cells")

# ── per-sample counts ─────────────────────────────────────────────────────────
if SAMPLE_COL not in sub.obs.columns:
    raise ValueError(f"Column '{SAMPLE_COL}' not in adata.obs.")

obs = sub.obs[[SAMPLE_COL]].copy()
for g, m in gene_masks.items():
    obs[g] = m.astype(int)
obs["Dcx+Prox1+Calb1"] = all_calb1.astype(int)
obs["Dcx+Prox1+Calb2"] = all_calb2.astype(int)
obs["total"] = 1

samples_present = [s for s in SAMPLE_ORDER if s in obs[SAMPLE_COL].values]
counts = obs.groupby(SAMPLE_COL).sum().reindex(samples_present).fillna(0).astype(int)

print("\nPer-sample counts:")
print(counts.to_string())

# ── figure: grouped bar chart per sample ──────────────────────────────────────
categories = GENES + ["Dcx+Prox1+Calb1", "Dcx+Prox1+Calb2"]
n_samples  = len(samples_present)
n_cats     = len(categories)

x        = np.arange(n_samples)
bar_w    = 0.13
offsets  = np.linspace(-(n_cats - 1) / 2 * bar_w, (n_cats - 1) / 2 * bar_w, n_cats)

fig, axes = plt.subplots(
    2, 1,
    figsize=(max(10, n_samples * 1.5), 7),
    gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
)
ax_bar = axes[0]
ax_pct = axes[1]

for i, cat in enumerate(categories):
    vals = [counts.loc[s, cat] if s in counts.index else 0 for s in samples_present]
    totals = [counts.loc[s, "total"] if s in counts.index else 1 for s in samples_present]
    ax_bar.bar(x + offsets[i], vals, width=bar_w, color=COLORS[cat], label=cat, zorder=3)

    pcts = [100 * v / t if t > 0 else 0 for v, t in zip(vals, totals)]
    ax_pct.bar(x + offsets[i], pcts, width=bar_w, color=COLORS[cat], zorder=3)

ax_bar.set_title(
    f"Marker expression in '{found}' cells\n"
    f"(Dcx, Prox1, Calb1, Calb2 — individually and co-expressed with Dcx+Prox1)",
    fontsize=11, fontweight="bold", pad=8,
)
ax_bar.set_ylabel("Cell count", fontsize=9)
ax_bar.set_xticks(x)
ax_bar.set_xticklabels([])
ax_bar.tick_params(axis="x", bottom=False)
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)
ax_bar.legend(fontsize=8.5, frameon=False, ncol=2)
ax_bar.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)

# annotate total cells per sample above bars
for xi, s in enumerate(samples_present):
    n = counts.loc[s, "total"] if s in counts.index else 0
    ax_bar.text(xi, ax_bar.get_ylim()[1] * 0.01, f"n={n}",
                ha="center", va="bottom", fontsize=7.5, color="#555555")

ax_pct.set_ylabel("% of cells", fontsize=9)
ax_pct.set_xticks(x)
ax_pct.set_xticklabels(samples_present, fontsize=9)
ax_pct.spines["top"].set_visible(False)
ax_pct.spines["right"].set_visible(False)
ax_pct.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)

COLOR_OIL  = "#E67E22"
COLOR_CORT = "#27AE60"
for tick_label, s in zip(ax_pct.get_xticklabels(), samples_present):
    if s in CORT_SAMPLES:
        tick_label.set_color(COLOR_CORT)
        tick_label.set_fontweight("bold")
    elif s in OIL_SAMPLES:
        tick_label.set_color(COLOR_OIL)
        tick_label.set_fontweight("bold")

# dashed separator between oil / cort groups
n_oil = sum(1 for s in samples_present if s in OIL_SAMPLES)
if 0 < n_oil < n_samples:
    for ax in [ax_bar, ax_pct]:
        ax.axvline(n_oil - 0.5, color="#CCCCCC", linewidth=1, linestyle="--", zorder=0)
    ax_bar.annotate("Oil",  xy=(n_oil / 2 - 0.5, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(0, 6), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8.5, color=COLOR_OIL, fontweight="bold")
    n_cort = n_samples - n_oil
    ax_bar.annotate("Cort", xy=(n_oil + n_cort / 2 - 0.5, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(0, 6), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8.5, color=COLOR_CORT, fontweight="bold")

os.makedirs(OUT_DIR, exist_ok=True)
out_path = os.path.join(OUT_DIR, OUT_FILE)
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nSaved: {out_path}")
plt.close()
