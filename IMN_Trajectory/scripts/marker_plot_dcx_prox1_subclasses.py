"""
Marker plot of the most differentially expressed genes per subclass,
restricted to cells that co-express Dcx and Prox1.

Produces two plots:
  1. marker_plot_dcx_prox1_subclasses_all.png
       — all subclasses retained regardless of sample representation.
  2. marker_plot_dcx_prox1_subclasses_filtered.png
       — only subclasses present in all 3 replicates of every condition are retained.

Also writes:
  out/subclass_replicate_representation.csv
       — per-subclass, per-condition replicate counts; flags condition-specific
         and under-represented subclasses.

Gene lookup follows the same logic as interactive_umap.py:
  - exact case-insensitive match on adata.var[sym_col], then partial fallback.

Usage:
    conda run -n anndata_env python scripts/marker_plot_dcx_prox1_subclasses.py
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")

# ── config ─────────────────────────────────────────────────────────────────────
DATA         = "data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
OUT_DIR      = "out"
os.makedirs(OUT_DIR, exist_ok=True)

FILTER_GENES   = ["Dcx", "Prox1"]   # cells must be positive for both
SUBCLASS_COL   = "subclass_name"
SAMPLE_COL     = "sample"
CONDITION_COL  = "treatment"        # column distinguishing experimental conditions
N_REPLICATES   = 3                  # required replicates per condition to keep a subclass
N_TOP_GENES    = 5                  # top marker genes per subclass to display

# ── load ───────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad(DATA)
print(f"  Shape: {adata.shape}")

# ── locate gene_symbol column (same logic as interactive_umap.py) ──────────────
sym_cols = [c for c in adata.var.columns if "gene_symbol" in c.lower()]
if not sym_cols:
    raise ValueError(f"No gene_symbol column found. Available: {adata.var.columns.tolist()}")
sym_col  = sym_cols[0]
print(f"  Using var column: '{sym_col}'")
symbols  = adata.var[sym_col].astype(str)


def resolve_gene(name):
    """Return (var_name, display_label) using the same lookup as interactive_umap.py."""
    mask = symbols.str.lower() == name.lower()
    hits = adata.var_names[mask].tolist()
    if not hits:
        mask = symbols.str.lower().str.contains(name.lower(), na=False)
        hits = adata.var_names[mask].tolist()
    if not hits:
        raise ValueError(f"Gene '{name}' not found in adata.var['{sym_col}'].")
    if len(hits) > 1:
        print(f"  WARNING: '{name}' matched {len(hits)} entries: {hits} — using first.")
    label = symbols[hits[0]]
    print(f"  Resolved '{name}' → var_name: {hits[0]}  symbol: {label}")
    return hits[0], label


def get_expr(var_name):
    idx = list(adata.var_names).index(var_name)
    col = adata.X[:, idx]
    return np.asarray(col.todense()).flatten() if sp.issparse(col) else np.asarray(col).flatten()


# ── build Dcx+ / Prox1+ cell mask ─────────────────────────────────────────────
print("\nFiltering cells...")
cell_mask = np.ones(adata.n_obs, dtype=bool)
for gene in FILTER_GENES:
    var_name, label = resolve_gene(gene)
    expr      = get_expr(var_name)
    cell_mask &= expr > 0
    print(f"  {label}+: {(expr > 0).sum()} cells")

n_pos = cell_mask.sum()
print(f"\n  Dcx+ & Prox1+ cells: {n_pos} / {adata.n_obs} ({n_pos/adata.n_obs*100:.2f}%)")

# ── base subset (all subclasses) ───────────────────────────────────────────────
sub_all = adata[cell_mask].copy()
sub_all.obs[SUBCLASS_COL] = sub_all.obs[SUBCLASS_COL].astype("category")
print(f"  Subset shape (all subclasses): {sub_all.shape}")


def make_dotplot(sub, label, out_filename):
    """Run DE and save a dot plot for the given AnnData subset."""
    n_groups = sub.obs[SUBCLASS_COL].nunique()
    print(f"\n[{label}] {n_groups} subclass(es), {sub.n_obs} cells")
    print(f"  Subclasses:")
    for sc_name in sub.obs[SUBCLASS_COL].cat.categories:
        n = (sub.obs[SUBCLASS_COL] == sc_name).sum()
        print(f"    {sc_name:<40} {n:>5}")

    # rename var_names to gene symbols for readable axis labels
    s = sub.copy()
    s.var["original_var_name"] = s.var_names.tolist()
    s.var_names                = s.var[sym_col].astype(str).tolist()
    s.var_names_make_unique()

    print(f"  Running rank_genes_groups (Wilcoxon)...")
    sc.tl.rank_genes_groups(
        s,
        groupby=SUBCLASS_COL,
        method="wilcoxon",
        n_genes=N_TOP_GENES * 2,
        use_raw=False,
        key_added="rank_genes",
    )

    print(f"  Plotting...")
    import matplotlib.pyplot as plt
    sc.pl.rank_genes_groups_dotplot(
        s,
        groupby=SUBCLASS_COL,
        n_genes=N_TOP_GENES,
        key="rank_genes",
        show=False,
        standard_scale="var",
        colorbar_title="Mean expression\n(scaled per gene)",
        size_title="Fraction of cells",
        figsize=(max(10, N_TOP_GENES * n_groups * 0.35), max(5, n_groups * 0.6 + 2)),
    )

    fig = plt.gcf()
    fig.suptitle(
        f"Top {N_TOP_GENES} marker genes per subclass — {label}\n"
        f"(Dcx+ / Prox1+ cells, Wilcoxon, n={sub.n_obs})",
        fontsize=11,
        y=1.02,
    )

    out_path = os.path.join(OUT_DIR, out_filename)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ── filter: keep only subclasses with ≥N_REPLICATES samples per condition ──────
print(f"\n  Building replicate-representation table (condition col: '{CONDITION_COL}')...")

obs = sub_all.obs[[SUBCLASS_COL, CONDITION_COL, SAMPLE_COL]].copy()
conditions     = obs[CONDITION_COL].unique().tolist()
all_subclasses = obs[SUBCLASS_COL].cat.categories.tolist()

# count unique samples per subclass × condition
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

csv_path = os.path.join(OUT_DIR, "subclass_replicate_representation_dcx_prox1.csv")
rep_counts.to_csv(csv_path)
print(f"  Saved representation table: {csv_path}")

keep    = rep_counts[rep_counts["status"] == "keep"].index.tolist()
dropped = rep_counts[rep_counts["status"] != "keep"].index.tolist()
print(f"\n  Subclasses kept  (≥{N_REPLICATES} replicates in every condition): {len(keep)}")
print(f"  Subclasses dropped: {len(dropped)}")
for sc_name in dropped:
    print(f"    {sc_name:<50}  [{rep_counts.loc[sc_name, 'status']}]  "
          + "  ".join(f"{c}={rep_counts.loc[sc_name, c]}" for c in conditions))

sub_filt = sub_all[sub_all.obs[SUBCLASS_COL].isin(keep)].copy()
sub_filt.obs[SUBCLASS_COL] = sub_filt.obs[SUBCLASS_COL].astype("category")

# ── version 1: all subclasses (filtered to ≥N_REPLICATES per condition) ────────
make_dotplot(
    sub_filt,
    label=f"all retained subclasses (≥{N_REPLICATES} replicates/condition)",
    out_filename="marker_plot_dcx_prox1_subclasses_all.png",
)

# ── version 2: subclasses present in ALL replicates of every condition ──────────
n_samples = sub_all.obs[SAMPLE_COL].nunique()
samples_per_subclass = (
    sub_filt.obs.groupby(SUBCLASS_COL)[SAMPLE_COL].nunique()
)
keep_all = samples_per_subclass[samples_per_subclass == n_samples].index.tolist()
sub_all_samples = sub_filt[sub_filt.obs[SUBCLASS_COL].isin(keep_all)].copy()
sub_all_samples.obs[SUBCLASS_COL] = sub_all_samples.obs[SUBCLASS_COL].astype("category")

make_dotplot(
    sub_all_samples,
    label=f"present in all {n_samples} samples",
    out_filename="marker_plot_dcx_prox1_subclasses_filtered.png",
)
