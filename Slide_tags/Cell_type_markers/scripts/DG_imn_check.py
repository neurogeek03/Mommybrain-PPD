import matplotlib
matplotlib.use('Agg')

from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path(__file__).resolve().parents[1]
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# =================== PARAMS ===================
group_col      = "subclass_name"
GENE_SYMBOL_COL = 'gene_symbol'

TARGET_TYPES = [
    '037 DG Glut',
    '038 DG-PIR Ex IMN',
    '045 OB-STR-CTX Inh IMN',
]

RANK_POOL_SIZE = 50   # Wilcoxon candidate pool per cell type
N_SHARED       = 10   # shared markers to show
N_UNIQUE       = 10   # unique markers per cell type to show

# =================== INPUT ===================
ad_file = '/scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'

adata = sc.read_h5ad(ad_file)
print(adata)

# =================== FILTER TO TARGET TYPES ===================
obs_unique = sorted(adata.obs[group_col].unique())
print(f"\nAll unique '{group_col}' values in file ({len(obs_unique)}):")
for v in obs_unique:
    print(f"  {v}")

missing = [ct for ct in TARGET_TYPES if ct not in obs_unique]
if missing:
    raise ValueError(f"Target cell types not found in '{group_col}': {missing}")

adata = adata[adata.obs[group_col].isin(TARGET_TYPES)].copy()
print(f"\nFiltered to {adata.n_obs} cells across {adata.obs[group_col].nunique()} cell types")
print(adata.obs[group_col].value_counts().to_string())

# =================== GENE SYMBOLS ===================
if GENE_SYMBOL_COL in adata.var.columns:
    adata.var_names = adata.var[GENE_SYMBOL_COL].astype(str)
    adata.var_names_make_unique()
    print(f"\nvar_names set to: {GENE_SYMBOL_COL}")
else:
    print(f"\nWARNING: '{GENE_SYMBOL_COL}' not found in adata.var — keeping existing var_names")
    print(f"  Available columns: {adata.var.columns.tolist()}")

# =================== PREPROCESSING ===================
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
print("Normalize + log1p done.")

# Fix categorical order
adata.obs[group_col] = pd.Categorical(
    adata.obs[group_col], categories=TARGET_TYPES, ordered=True
)

# =================== RANK GENES (Wilcoxon) ===================
print("\nRunning rank_genes_groups (Wilcoxon)...")
sc.tl.rank_genes_groups(adata, groupby=group_col, method='wilcoxon', use_raw=False)

result = adata.uns['rank_genes_groups']
candidates = {}  # {cell_type: [gene, ...]}
scores     = {}  # {cell_type: [score, ...]}

for ct in TARGET_TYPES:
    candidates[ct] = result['names'][ct][:RANK_POOL_SIZE].tolist()
    scores[ct]     = result['scores'][ct][:RANK_POOL_SIZE].tolist()
    print(f"  {ct}: top 5 candidates = {candidates[ct][:5]}")

# =================== SHARED vs UNIQUE ===================

# Shared: gene in top-N pool of ALL 3 cell types
shared_set = (
    set(candidates[TARGET_TYPES[0]]) &
    set(candidates[TARGET_TYPES[1]]) &
    set(candidates[TARGET_TYPES[2]])
)
print(f"\nGenes in top-{RANK_POOL_SIZE} of ALL 3 types: {len(shared_set)}")

# Rank shared genes by average Wilcoxon score across the 3 types
shared_avg_score = {}
for gene in shared_set:
    avg = np.mean([
        scores[ct][candidates[ct].index(gene)]
        for ct in TARGET_TYPES
    ])
    shared_avg_score[gene] = avg

shared_genes = sorted(shared_avg_score, key=shared_avg_score.get, reverse=True)[:N_SHARED]
print(f"Top {N_SHARED} shared genes: {shared_genes}")

# Unique: gene in top-N of EXACTLY 1 cell type (absent from both others)
unique_genes = {}
for ct in TARGET_TYPES:
    others = [c for c in TARGET_TYPES if c != ct]
    other_union = set(candidates[others[0]]) | set(candidates[others[1]])
    unique_to_ct = [g for g in candidates[ct] if g not in other_union][:N_UNIQUE]
    unique_genes[ct] = unique_to_ct
    print(f"Unique to '{ct}' (top {N_UNIQUE}): {unique_to_ct}")

# =================== GENE LIST ASSEMBLY ===================
# Order: shared block, then unique block per cell type
gene_groups = {'Shared': shared_genes}
for ct in TARGET_TYPES:
    short_label = ct.split(' ', 1)[1]  # drop numeric prefix for legend readability
    gene_groups[f"Unique\n{short_label}"] = unique_genes[ct]

all_genes = [g for genes in gene_groups.values() for g in genes]
print(f"\nTotal genes in dotplot: {len(all_genes)}")

# =================== DOTPLOT ===================
sc.settings.figdir = str(output_base)

dp = sc.pl.dotplot(
    adata,
    var_names=gene_groups,
    groupby=group_col,
    use_raw=False,
    standard_scale='var',   # scale each gene to [0,1] across the 3 groups
    swap_axes=True,         # genes on y-axis, cell types on x-axis
    return_fig=True,
    show=False,
)

pdf_path = output_base / f"DG_imn_check_{stamp}.pdf"
png_path = output_base / f"DG_imn_check_{stamp}.png"
dp.savefig(str(pdf_path), bbox_inches='tight')
dp.savefig(str(png_path), bbox_inches='tight', dpi=150)
print(f"\nSaved: {pdf_path.name}")
print(f"Saved: {png_path.name}")

# =================== SUMMARY CSV ===================
rows = []
for category, genes in gene_groups.items():
    for gene in genes:
        row = {'category': category, 'gene': gene}
        for ct in TARGET_TYPES:
            in_pool = gene in candidates[ct]
            rank    = candidates[ct].index(gene) + 1 if in_pool else None
            score   = scores[ct][candidates[ct].index(gene)] if in_pool else None
            row[f"rank_{ct}"]  = rank
            row[f"score_{ct}"] = score
        rows.append(row)

summary_df = pd.DataFrame(rows)
csv_path = output_base / f"DG_imn_check_markers_{stamp}.csv"
summary_df.to_csv(csv_path, index=False)
print(f"Saved: {csv_path.name}")

print(f"\n------ Script completed successfully ------")