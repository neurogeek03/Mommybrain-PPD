"""
correlate.py

Compute per-cell-type Pearson r between slide-seq and slide-tags pseudobulk
mean expression profiles, using shared genes (Ensembl IDs).

Usage:
  uv run scripts/correlate.py              # all shared genes
  uv run scripts/correlate.py --de-only    # DE genes only (FDR<0.05, |logFC|>1)

Outputs (in output/ or output/de_only/ when --de-only):
  correlation_per_celltype.csv
  summary_barplot.png
  per_celltype_gene_scatter.pdf
  final_scatter_celltype.png
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import anndata as ad

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--de-only", action="store_true",
                    help="Restrict to DE genes per cell type (FDR<0.05, |logFC|>1)")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / ("output/de_only" if args.de_only else "output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SLIDESEQ_H5AD = DATA_DIR / "All_RCTD_types_B03_B14_filtered.h5ad"
SLIDETAGS_H5AD = DATA_DIR / "PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
COLOR_CSV = Path("/scratch/mfafouti/CAN-2026-poster/data/cluster_annotation_term.csv")
EDGER_DIR = Path("/scratch/mfafouti/CAN-2026-poster/data/EdgeR/OIL_vs_CORT")

ASTRO_PREFIXES = {"318", "319", "321"}
ASTRO_COLOR = "#AD5CCC"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def prefix_from_slideseq(ct: str) -> str:
    """Extract 3-digit numeric prefix from a slide-seq cell_type string.
    e.g. '016_CA1_ProS_Glut' -> '016'
    'Astrocytes' -> 'Astrocytes'  (special case, already merged)
    """
    m = re.match(r"^(\d{3})", ct)
    return m.group(1) if m else ct


def normalize_slideseq_name(ct: str) -> str:
    """Return a canonical key for a slide-seq cell type.
    For 'Astrocytes' the key is 'Astrocytes'.
    For numeric types the key is the 3-digit prefix.
    """
    if ct == "Astrocytes":
        return "Astrocytes"
    return prefix_from_slideseq(ct)


def normalize_slidetags_name(subclass: str) -> str:
    """Return canonical key for a slide-tags subclass_name.
    subclass_name format: '016 CA1-ProS Glut'
    Returns '016' (3-digit prefix), or 'Astrocytes' for 318/319/321.
    """
    m = re.match(r"^(\d{3})", subclass)
    if m:
        pfx = m.group(1)
        if pfx in ASTRO_PREFIXES:
            return "Astrocytes"
        return pfx
    return subclass


def build_color_map(color_csv: Path) -> dict:
    """Build prefix -> hex color dict from the subclass rows of cluster_annotation_term.csv."""
    df = pd.read_csv(color_csv)
    subc = df[df["cluster_annotation_term_set_name"] == "subclass"].copy()
    color_map = {}
    for _, row in subc.iterrows():
        name = str(row["name"])
        m = re.match(r"^(\d{3})", name)
        if m:
            color_map[m.group(1)] = row["color_hex_triplet"]
    return color_map


def load_de_gene_sets(edger_dir: Path, symbol_to_ensembl: dict,
                      fdr_thresh: float = 0.05, lfc_thresh: float = 1.0) -> dict:
    """Return {prefix -> set of Ensembl IDs} for each EdgeR result file.
    Astrocytes key is 'Astrocytes' (union of 318 + 319 files).
    """
    de_sets = {}
    astro_union: set = set()
    for fpath in sorted(edger_dir.glob("*_edgeR_results.tsv")):
        m = re.match(r"^(\d{3})", fpath.name)
        if not m:
            continue
        prefix = m.group(1)
        df = pd.read_csv(fpath, sep="\t", index_col=0)
        sig = df[(df["FDR"] < fdr_thresh) & (df["logFC"].abs() > lfc_thresh)]
        ensembl_ids = {symbol_to_ensembl[g] for g in sig.index if g in symbol_to_ensembl}
        if prefix in ASTRO_PREFIXES:
            astro_union |= ensembl_ids
        else:
            de_sets[prefix] = ensembl_ids
    if astro_union:
        de_sets["Astrocytes"] = astro_union
    return de_sets


def pseudobulk_mean(adata: ad.AnnData, cell_mask: np.ndarray) -> np.ndarray:
    """Compute mean expression across selected cells/spots.

    Works for both dense and sparse X.
    Returns a 1-D array of length n_genes.
    """
    X_sub = adata.X[cell_mask]
    if hasattr(X_sub, "toarray"):
        X_sub = X_sub.toarray()
    return np.asarray(X_sub.mean(axis=0)).ravel()


# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
print("Loading slide-seq...")
ss = ad.read_h5ad(SLIDESEQ_H5AD)
print(f"  Slide-seq: {ss.n_obs} spots x {ss.n_vars} genes")

print("Loading slide-tags...")
st = ad.read_h5ad(SLIDETAGS_H5AD)
print(f"  Slide-tags: {st.n_obs} cells x {st.n_vars} genes")

# ---------------------------------------------------------------------------
# 2. Subset both to shared genes
# ---------------------------------------------------------------------------
shared_genes = ss.var_names.intersection(st.var_names)
print(f"Shared genes: {len(shared_genes)}")

ss = ss[:, shared_genes].copy()
st = st[:, shared_genes].copy()

# ---------------------------------------------------------------------------
# 2b. Build symbol -> Ensembl map and load DE gene sets (if --de-only)
# ---------------------------------------------------------------------------
symbol_to_ensembl = {
    sym: eid
    for eid, sym in st.var["gene_symbols"].dropna().items()
    if eid in shared_genes
}

if args.de_only:
    de_sets = load_de_gene_sets(EDGER_DIR, symbol_to_ensembl)
    print(f"Loaded DE gene sets for {len(de_sets)} cell types")
    for k, s in sorted(de_sets.items()):
        print(f"  {k}: {len(s)} DE genes (after Ensembl mapping)")
else:
    de_sets = {}

# ---------------------------------------------------------------------------
# 3. Create normalized cell_type column in slide-tags
# ---------------------------------------------------------------------------
st.obs["cell_type"] = st.obs["subclass_name"].apply(normalize_slidetags_name)

# ---------------------------------------------------------------------------
# 4. Find shared cell types
# ---------------------------------------------------------------------------
ss_types = {normalize_slideseq_name(ct): ct for ct in ss.obs["cell_type"].unique()}
st_types = set(st.obs["cell_type"].unique())

all_shared = sorted(set(ss_types.keys()) & st_types)

# When --de-only, restrict to cell types that have an EdgeR file
if args.de_only:
    shared_keys = [k for k in all_shared if k in de_sets]
    print(f"Shared cell types with EdgeR results: {len(shared_keys)} / {len(all_shared)}")
else:
    shared_keys = all_shared
    print(f"Shared cell types: {len(shared_keys)}")
for k in shared_keys:
    print(f"  {k}  ->  ss: {ss_types[k]}")

# ---------------------------------------------------------------------------
# 5-6. Pseudobulk mean + Pearson r
# ---------------------------------------------------------------------------
results = []
pseudobulks = {}  # key -> (mean_st, mean_ss) for scatter plots

for key in shared_keys:
    ss_label = ss_types[key]

    ss_mask = (ss.obs["cell_type"] == ss_label).values
    st_mask = (st.obs["cell_type"] == key).values

    n_ss = ss_mask.sum()
    n_st = st_mask.sum()

    mean_ss = pseudobulk_mean(ss, ss_mask)
    mean_st = pseudobulk_mean(st, st_mask)

    if args.de_only:
        de_ensembl = de_sets.get(key, set()) & set(shared_genes)
        gene_idx = np.array([i for i, g in enumerate(shared_genes) if g in de_ensembl])
        if len(gene_idx) < 2:
            print(f"  {ss_label}: skipped (only {len(gene_idx)} DE gene(s) after mapping)")
            continue
        mean_ss = mean_ss[gene_idx]
        mean_st = mean_st[gene_idx]

    r, pval = stats.pearsonr(mean_ss, mean_st)

    pseudobulks[key] = (mean_st, mean_ss, ss_label, n_ss, n_st)

    results.append({
        "cell_type": ss_label,
        "cell_type_key": key,
        "pearson_r": r,
        "pearson_pval": pval,
        "n_slideseq_spots": int(n_ss),
        "n_slidetags_cells": int(n_st),
        "mean_st": float(mean_st.mean()),
        "mean_ss": float(mean_ss.mean()),
    })
    print(f"  {ss_label}: r={r:.4f}  n_ss={n_ss}  n_st={n_st}")

# ---------------------------------------------------------------------------
# 7. Save CSV
# ---------------------------------------------------------------------------
results_df = pd.DataFrame(results)
out_csv = OUTPUT_DIR / "correlation_per_celltype.csv"
results_df[
    ["cell_type", "pearson_r", "pearson_pval", "n_slideseq_spots", "n_slidetags_cells"]
].to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")

# ---------------------------------------------------------------------------
# 8. Summary bar chart
# ---------------------------------------------------------------------------
color_map = build_color_map(COLOR_CSV)

# Sort by r descending
plot_df = results_df.sort_values("pearson_r", ascending=False).reset_index(drop=True)

# Assign colors
def get_color(row):
    key = row["cell_type_key"]
    if key == "Astrocytes":
        return ASTRO_COLOR
    return color_map.get(key, "#AAAAAA")

plot_df["color"] = plot_df.apply(get_color, axis=1)

n_types = len(plot_df)
fig_height = max(6, n_types * 0.38)
fig, ax = plt.subplots(figsize=(10, fig_height))

bars = ax.barh(
    y=range(n_types),
    width=plot_df["pearson_r"],
    color=plot_df["color"],
    edgecolor="white",
    linewidth=0.4,
    height=0.75,
)

# Labels on right end of each bar showing n_spots / n_cells
for i, row in plot_df.iterrows():
    r_val = row["pearson_r"]
    label = f"n_ss={row['n_slideseq_spots']:,}  n_st={row['n_slidetags_cells']:,}"
    x_pos = r_val + 0.003 if r_val >= 0 else r_val - 0.003
    ha = "left" if r_val >= 0 else "right"
    ax.text(
        x_pos,
        i,
        label,
        va="center",
        ha=ha,
        fontsize=5.5,
        color="#444444",
    )

ax.set_yticks(range(n_types))
ax.set_yticklabels(plot_df["cell_type"], fontsize=7)
ax.invert_yaxis()  # highest r at top

gene_subtitle = "DE genes only (FDR<0.05, |logFC|>1)" if args.de_only else "all shared genes"
ax.set_xlabel("Pearson r  (slide-seq vs. slide-tags pseudobulk)", fontsize=9)
ax.set_title(
    f"Per-cell-type transcriptomic correlation — {gene_subtitle}\n(Slide-seq B03+B14  vs  Slide-tags all samples)",
    fontsize=10,
    pad=10,
)

ax.axvline(0, color="black", linewidth=0.7, linestyle="--")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", length=0)

# Legend: grey = unmatched in color file
legend_handles = [
    mpatches.Patch(facecolor="#AD5CCC", label="Astrocytes"),
    mpatches.Patch(facecolor="#AAAAAA", label="No color annotation"),
]
ax.legend(handles=legend_handles, fontsize=7, loc="lower right")

plt.tight_layout()
out_png = OUTPUT_DIR / "summary_barplot.png"
plt.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_png}")

# ---------------------------------------------------------------------------
# 9. Per-cell-type scatter plots (one gene per point, x=slide-tags, y=slide-seq)
#    Saved as a multi-page PDF, 12 panels per page (3 cols x 4 rows)
# ---------------------------------------------------------------------------
from matplotlib.backends.backend_pdf import PdfPages

COLS = 3
ROWS = 4
PER_PAGE = COLS * ROWS

sorted_keys = plot_df["cell_type_key"].tolist()  # already sorted by r desc

out_pdf = OUTPUT_DIR / "per_celltype_gene_scatter.pdf"
with PdfPages(out_pdf) as pdf:
    page_keys = []
    for i, key in enumerate(sorted_keys):
        page_keys.append(key)
        if len(page_keys) == PER_PAGE or i == len(sorted_keys) - 1:
            fig, axes = plt.subplots(ROWS, COLS, figsize=(14, 18))
            axes_flat = axes.ravel()
            for ax_idx, k in enumerate(page_keys):
                ax = axes_flat[ax_idx]
                mean_st, mean_ss, ss_label, n_ss, n_st = pseudobulks[k]
                row_data = results_df[results_df["cell_type_key"] == k].iloc[0]
                r_val = row_data["pearson_r"]
                color = ASTRO_COLOR if k == "Astrocytes" else color_map.get(k, "#AAAAAA")

                ax.scatter(mean_st, mean_ss, s=0.8, alpha=0.4, color=color, linewidths=0)
                lim_min = min(mean_st.min(), mean_ss.min())
                lim_max = max(mean_st.max(), mean_ss.max())
                ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=0.6, alpha=0.5)
                short_label = ss_label if len(ss_label) <= 28 else ss_label[:25] + "..."
                ax.set_title(f"{short_label}\nr={r_val:.3f}  n_ss={n_ss}  n_st={n_st}",
                             fontsize=5.5, pad=3)
                ax.set_xlabel("Slide-tags mean", fontsize=5)
                ax.set_ylabel("Slide-seq mean", fontsize=5)
                ax.tick_params(labelsize=4)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)

            # hide unused panels on last page
            for ax_idx in range(len(page_keys), PER_PAGE):
                axes_flat[ax_idx].set_visible(False)

            plt.suptitle(f"Per-cell-type: slide-tags vs slide-seq mean expression per gene ({gene_subtitle})",
                         fontsize=9, y=1.01)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            page_keys = []

print(f"Saved: {out_pdf}")

# ---------------------------------------------------------------------------
# 10. Final summary scatter: 1 point per cell type
#     x = mean expression across all genes in slide-tags
#     y = mean expression across all genes in slide-seq
# ---------------------------------------------------------------------------
final_df = results_df.copy()
final_df["color"] = final_df.apply(get_color, axis=1)

fig, ax = plt.subplots(figsize=(9, 8))

for _, row in final_df.iterrows():
    ax.scatter(row["mean_st"], row["mean_ss"], color=row["color"],
               s=60, edgecolors="white", linewidths=0.4, zorder=3)
    ax.annotate(
        row["cell_type_key"] if row["cell_type_key"] != "Astrocytes" else "Astrocytes",
        (row["mean_st"], row["mean_ss"]),
        fontsize=3.5, ha="left", va="bottom",
        xytext=(2, 2), textcoords="offset points", color="#333333",
    )

all_x = final_df["mean_st"]
all_y = final_df["mean_ss"]
lim_min = min(all_x.min(), all_y.min()) * 0.95
lim_max = max(all_x.max(), all_y.max()) * 1.05
ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=0.8, alpha=0.4, label="y = x")

ax.set_xlabel("Mean log-norm expression — Slide-tags", fontsize=10)
ax.set_ylabel("Mean log-norm expression — Slide-seq", fontsize=10)
ax.set_title(
    f"Mean expression per cell type: Slide-seq vs Slide-tags\n({gene_subtitle}, each point = one cell type)",
    fontsize=10, pad=8,
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(fontsize=8)

plt.tight_layout()
out_final = OUTPUT_DIR / "final_scatter_celltype.png"
plt.savefig(out_final, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_final}")
