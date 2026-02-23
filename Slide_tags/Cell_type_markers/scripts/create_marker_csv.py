from pathlib import Path
import re
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path(__file__).resolve().parents[1]
print(f"project path: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# =================== PARAMS ===================
group_col = "subclass_name"

RANK_POOL_SIZE = 50             # candidate genes per subclass from rank_genes_groups
TOP_N_GENES    = 2              # final marker genes per subclass (after specificity re-ranking)
MIN_SPEC_SCORE = 0.0            # specificity floor (0 = no filter)
GENE_SYMBOL_COL = 'gene_symbol' # adata.var column with human-readable symbols

# --- Subclass filter (Step 2): read from shared subclasses CSV ---
SUBCLASS_CSV = project_path / 'subclasses' / 'shared_subclasses.csv'
if SUBCLASS_CSV.exists():
    _sub_df = pd.read_csv(SUBCLASS_CSV)
    SELECTED_SUBCLASSES = _sub_df['subclass_name'].tolist()
    print(f"Loaded {len(SELECTED_SUBCLASSES)} subclasses from {SUBCLASS_CSV.name}")
else:
    SELECTED_SUBCLASSES = []
    print(f"WARNING: {SUBCLASS_CSV} not found — using all subclasses")

# --- Curated marker genes (Step 4): fill after reviewing candidate genes ---
CURATED_MARKERS = {
    # Broad class markers
    # '__Glutamatergic': ['GeneA'],
    # '__GABAergic':     ['GeneB'],
    # '__Non-neuronal':  ['GeneC'],
    # Per-subclass markers (one each)
    # 'SubclassX_Glut':  ['GeneD'],
    # 'SubclassY_Gaba':  ['GeneE'],
    # ...
}

# =================== INPUT ===================
ad_file ='/scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/DE_after_mt_filter_108123_mincells_10_in_2_samples_slide_tags.h5ad'
# ad_file = '/scratch/mfafouti/Mommybrain/Slide_tags/Cell_type_markers/subclasses/data/filtered_combined_20260222_173834.h5ad'

# --- Read data ---
adata = sc.read_h5ad(ad_file)
print(adata)

# --- Diagnostic: unique values in group column ---
obs_unique = sorted(adata.obs[group_col].unique())
print(f"\nUnique values in '{group_col}': {len(obs_unique)}")
for v in obs_unique:
    print(f"  {v}")

# --- Replace Ensembl IDs in var_names with gene symbols ---
if GENE_SYMBOL_COL in adata.var.columns:
    adata.var_names = adata.var[GENE_SYMBOL_COL].astype(str)
    adata.var_names_make_unique()
    print(f"adata.var_names now using: {GENE_SYMBOL_COL}")
else:
    print(f"WARNING: '{GENE_SYMBOL_COL}' not found in adata.var.columns: {adata.var.columns.tolist()}")
    print("Keeping existing var_names (update GENE_SYMBOL_COL if needed)")

# =================== PREPROCESSING ===================
print("\nNormalizing and log-transforming data...")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
print("Preprocessing complete.")

# =================== VAR INSPECTION ===================
print("\n--- Inspecting adata.var ---")
print(f"adata.var columns: {adata.var.columns.tolist()}")
print(f"\nFirst 5 rows of adata.var:")
print(adata.var.head())
print(f"\nadata.var_names (index, first 10):")
print(adata.var_names[:10].tolist())
# =====================================================

groups = adata.obs[group_col].unique()

# =================== GROUP DIAGNOSTICS ===================
n_groups = adata.obs[group_col].nunique()
print(f"\nNumber of unique {group_col}: {n_groups}")
print(f"Groups included: {sorted(groups)}")

# =================== SUBCLASS EXPLORATION (Step 1) ===================

def get_broad_class(subclass_name):
    if subclass_name.endswith('Glut'):
        return 'Glutamatergic'
    elif subclass_name.endswith('Gaba'):
        return 'GABAergic'
    elif subclass_name.endswith('NN'):
        return 'Non-neuronal'
    return 'Unknown'

cell_counts = adata.obs[group_col].value_counts()
count_df = cell_counts.reset_index()
count_df.columns = [group_col, 'n_cells']
count_df['broad_class'] = count_df[group_col].apply(get_broad_class)
count_df = count_df.sort_values('n_cells', ascending=False)

print("\n--- Cell counts per broader class ---")
print(count_df.groupby('broad_class')['n_cells'].agg(['sum', 'count']))
print(f"\n--- All subclasses ranked by cell count (n={len(count_df)}) ---")
print(count_df.to_string(index=False))

count_df.to_csv(output_base / f"subclass_cell_counts_{stamp}.csv", index=False)
print(f"Saved subclass_cell_counts_{stamp}.csv")

# =================== SUBCLASS FILTER (Step 2) ===================

if SELECTED_SUBCLASSES:
    # Normalise spaces/underscores so CSV names match obs labels
    obs_labels = set(adata.obs[group_col].unique())
    if not obs_labels & set(SELECTED_SUBCLASSES):
        SELECTED_SUBCLASSES = [re.sub(r"[ /-]", "_", s) for s in SELECTED_SUBCLASSES]
        if not obs_labels & set(SELECTED_SUBCLASSES):
            SELECTED_SUBCLASSES = [s.replace('_', ' ') for s in SELECTED_SUBCLASSES]
    adata = adata[adata.obs[group_col].isin(SELECTED_SUBCLASSES)].copy()
    groups = adata.obs[group_col].unique()
    print(f"\nFiltered to {len(groups)} subclasses ({adata.n_obs} cells)")

# --- Define subclass display order (groups broad classes together) ---
broad_class_order = ['Glutamatergic', 'GABAergic', 'Non-neuronal', 'Unknown']
subclass_order = sorted(
    groups,
    key=lambda s: (broad_class_order.index(get_broad_class(s)), s)
)
adata.obs[group_col] = pd.Categorical(
    adata.obs[group_col], categories=subclass_order, ordered=True
)
print(f"\nSubclass display order ({len(subclass_order)}):")
for s in subclass_order:
    print(f"  {get_broad_class(s):16s}  {s}")

# =================== DIAGNOSTICS ===================

def check_log_normalized(adata):
    """
    Check if adata.X is log-normalized.

    Checks three independent signals:
    1. 'log1p' in adata.uns — written by sc.pp.log1p()
    2. max(adata.X) < 15 — log-normalized data is bounded; raw counts are not
    3. adata.raw is not None — informational (doesn't confirm X is log-norm)

    Returns bool, does not auto-normalize.
    """
    signal1 = 'log1p' in adata.uns

    # Get max value safely for sparse or dense
    if sp.issparse(adata.X):
        max_val = adata.X.max()
    else:
        max_val = np.max(adata.X) if adata.X.size > 0 else 0
    signal2 = max_val < 15

    signal3 = adata.raw is not None

    print("\n--- Log-Normalization Diagnostic ---")
    print(f"  'log1p' in adata.uns: {signal1}")
    print(f"  max(adata.X) < 15: {signal2} (max={max_val:.2f})")
    print(f"  adata.raw is not None: {signal3}")
    print(f"  Likely log-normalized: {signal1 or signal2}")
    print()

    return signal1 or signal2

is_log_norm = check_log_normalized(adata)
if not is_log_norm:
    print("WARNING: adata.X may not be log-normalized. sc.pl.dotplot requires log-norm data.")


# =================== RANK GENES ===================

print("\nRanking genes by group using scanpy...")
sc.tl.rank_genes_groups(
    adata,
    groupby=group_col,
    method='wilcoxon',
    use_raw=False,
)

# Extract candidate pool per group
print(f"\nExtracting top {RANK_POOL_SIZE} candidate genes per subclass...")
result = adata.uns['rank_genes_groups']
candidate_genes_dict = {}

for group in groups:
    top_genes = result['names'][group][:RANK_POOL_SIZE]
    candidate_genes_dict[group] = top_genes.tolist()

# Unique candidate genes across all subclasses
all_candidates = list({g for genes in candidate_genes_dict.values() for g in genes})
print(f"Total unique candidate genes across all subclasses: {len(all_candidates)}")

# =================== SPECIFICITY SCORING ===================

print("\nComputing specificity scores for candidate genes...")

# Subset to candidate genes only (memory-safe)
adata_cand = adata[:, all_candidates]

# Percent expressing per group (binary: expressed > 0)
pct_expr = pd.DataFrame(index=all_candidates, columns=list(groups), dtype=float)
for group in groups:
    mask = adata_cand.obs[group_col] == group
    n_cells = mask.sum()
    if n_cells == 0:
        pct_expr[group] = 0.0
        continue
    group_X = adata_cand[mask].X
    if sp.issparse(group_X):
        pct_expr[group] = np.asarray((group_X > 0).sum(axis=0)).flatten() / n_cells
    else:
        pct_expr[group] = (group_X > 0).sum(axis=0) / n_cells

# Specificity = pct_expr_in_group - max(pct_expr_in_all_other_groups)
specificity = pd.DataFrame(index=all_candidates, columns=list(groups), dtype=float)
for group in groups:
    other_max = pct_expr.drop(columns=group).max(axis=1)
    specificity[group] = pct_expr[group] - other_max

# Save full specificity table
spec_csv = output_base / f"candidate_specificity_scores_{stamp}.csv"
specificity.to_csv(spec_csv)
print(f"Saved specificity scores: {spec_csv.name}")

# =================== TOP MARKERS BY SPECIFICITY ===================

# For each subclass, re-rank its candidate pool by specificity and pick top N
marker_genes_dict = {}
print(f"\nSelecting top {TOP_N_GENES} genes per subclass by specificity:")
for group in subclass_order:
    candidates = candidate_genes_dict.get(group, [])
    group_spec = specificity.loc[candidates, group].sort_values(ascending=False)
    if MIN_SPEC_SCORE > 0:
        group_spec = group_spec[group_spec >= MIN_SPEC_SCORE]
    top = group_spec.head(TOP_N_GENES)
    marker_genes_dict[group] = top.index.tolist()
    broad = get_broad_class(group)
    print(f"  [{broad}] {group}: {list(zip(top.index, top.values.round(3)))}")

# =================== GENE LIST ASSEMBLY ===================

if CURATED_MARKERS:
    marker_genes_list = []
    for genes in CURATED_MARKERS.values():
        marker_genes_list.extend(genes)
    marker_genes_list = list(dict.fromkeys(marker_genes_list))
    print(f"\nUsing CURATED_MARKERS: {len(marker_genes_list)} unique genes")
else:
    # Flatten in subclass display order for diagonal pattern
    marker_genes_list = []
    for group in subclass_order:
        for gene in marker_genes_dict.get(group, []):
            if gene not in marker_genes_list:
                marker_genes_list.append(gene)
    print(f"\nUsing specificity-ranked markers (diagonal order): {len(marker_genes_list)} unique genes")


# =================== DOTPLOT DATA EXPORT ===================

print("\nComputing per-(gene, subclass) summary for dotplot data...")

adata_markers = adata[:, marker_genes_list]

rows = []
for i, group in enumerate(subclass_order):
    mask = adata_markers.obs[group_col] == group
    n_cells = mask.sum()
    if n_cells == 0:
        for gene in marker_genes_list:
            rows.append({
                'subclass': group,
                'gene': gene,
                'mean_expr': 0.0,
                'pct_expr': 0.0,
                'specificity': 0.0,
                'broad_class': get_broad_class(group),
                'subclass_order': i,
            })
        continue

    group_X = adata_markers[mask].X
    for j, gene in enumerate(marker_genes_list):
        col = group_X[:, j]
        if sp.issparse(col):
            col = np.asarray(col.todense()).flatten()
        else:
            col = np.asarray(col).flatten()

        mean_val = float(col.mean())
        pct_val = float((col > 0).sum() / n_cells)
        spec_val = float(specificity.at[gene, group]) if gene in specificity.index else 0.0

        rows.append({
            'subclass': group,
            'gene': gene,
            'mean_expr': mean_val,
            'pct_expr': pct_val,
            'specificity': spec_val,
            'broad_class': get_broad_class(group),
            'subclass_order': i,
        })

dotplot_df = pd.DataFrame(rows)
dotplot_csv = output_base / f"marker_dotplot_data.csv"
dotplot_df.to_csv(dotplot_csv, index=False)
print(f"Saved dotplot data: {dotplot_csv}")
print(f"  Shape: {dotplot_df.shape[0]} rows  ({len(marker_genes_list)} genes x {len(subclass_order)} subclasses)")

print(f"\n------ Script completed successfully ------")

