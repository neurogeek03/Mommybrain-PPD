import matplotlib
matplotlib.use('Agg')

from pathlib import Path
import scanpy as sc
import pandas as pd
from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path(__file__).resolve().parents[1]
output_base  = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)

# =================== PARAMS ===================
group_col       = "subclass_name"
GENE_SYMBOL_COL = 'gene_symbol'

TARGET_TYPES = [
    '037 DG Glut',
    '038 DG-PIR Ex IMN',
    '045 OB-STR-CTX Inh IMN',
]

# Curated marker genes per cell type
# NOTE: verify spellings — possible typos flagged with '?'
CURATED_MARKERS = {
    '037 DG Glut':            ['Prox1', 'Itoka', 'C1q13'],  # Itoka → Itpka? C1q13 → C1ql3?
    '038 DG-PIR Ex IMN':      ['Mex3a', 'Neurod1'],
    '045 OB-STR-CTX Inh IMN': ['DIx1', 'Zeb2'],             # DIx1 → Dlx1?
}

# =================== INPUT ===================
ad_file = '/scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad'

adata = sc.read_h5ad(ad_file)
print(adata)

# =================== FILTER TO TARGET TYPES ===================
obs_unique = sorted(adata.obs[group_col].unique())
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
    print(f"\nWARNING: '{GENE_SYMBOL_COL}' not found — keeping existing var_names")
    print(f"  Available columns: {adata.var.columns.tolist()}")

# =================== PREPROCESSING ===================
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
print("Normalize + log1p done.")

adata.obs[group_col] = pd.Categorical(
    adata.obs[group_col], categories=TARGET_TYPES, ordered=True
)

# =================== DCX CHECK ===================
DCX_GENE = 'Dcx'  # adjust capitalisation if needed (e.g. 'DCX')

print(f"\n--- DCX expression check ({DCX_GENE}) ---")
if DCX_GENE not in adata.var_names:
    # Case-insensitive fallback
    matches = [g for g in adata.var_names if g.lower() == DCX_GENE.lower()]
    if matches:
        DCX_GENE = matches[0]
        print(f"  Found as: {DCX_GENE}")
    else:
        print(f"  WARNING: '{DCX_GENE}' not found in var_names — skipping DCX check")
        DCX_GENE = None

if DCX_GENE:
    import scipy.sparse as sp
    dcx_vec = adata[:, DCX_GENE].X
    if sp.issparse(dcx_vec):
        dcx_vec = dcx_vec.toarray().flatten()
    else:
        dcx_vec = dcx_vec.flatten()

    dcx_series = pd.Series(dcx_vec, index=adata.obs_names)
    for ct in TARGET_TYPES:
        mask   = adata.obs[group_col] == ct
        vals   = dcx_series[mask]
        pct    = (vals > 0).mean() * 100
        mean_e = vals[vals > 0].mean() if (vals > 0).any() else 0.0
        print(f"  {ct}")
        print(f"    % cells expressing: {pct:.1f}%  ({(vals > 0).sum()}/{mask.sum()} cells)")
        print(f"    mean expr (expressing cells only): {mean_e:.3f}")

# =================== GENE VALIDATION ===================
all_curated = [g for genes in CURATED_MARKERS.values() for g in genes]
missing_genes = [g for g in all_curated if g not in adata.var_names]
found_genes   = [g for g in all_curated if g     in adata.var_names]

print(f"\nCurated genes found in data ({len(found_genes)}/{len(all_curated)}): {found_genes}")
if missing_genes:
    print(f"WARNING: genes not found in var_names — check spelling: {missing_genes}")

# Build gene_groups using only genes that exist in the data
gene_groups = {}
for ct, genes in CURATED_MARKERS.items():
    valid = [g for g in genes if g in adata.var_names]
    if valid:
        short_label = ct.split(' ', 1)[1]  # drop numeric prefix
        gene_groups[short_label] = valid

if DCX_GENE:
    gene_groups['DCX'] = [DCX_GENE]

if not gene_groups:
    raise RuntimeError("No curated genes found in data. Check gene names and GENE_SYMBOL_COL.")

all_plot_genes = [g for genes in gene_groups.values() for g in genes]
print(f"\nGene groups for dotplot:")
for label, genes in gene_groups.items():
    print(f"  {label}: {genes}")

# =================== DOTPLOT ===================
import matplotlib.pyplot as plt

dp = sc.pl.dotplot(
    adata,
    var_names=gene_groups,
    groupby=group_col,
    use_raw=False,
    standard_scale='var',   # scale each gene to [0,1] across groups
    swap_axes=True,         # genes on y-axis, cell types on x-axis
    return_fig=True,
    show=False,
)

# Rotate cell type labels on x-axis so they are visible
main_ax = dp.get_axes()['mainplot_ax']
main_ax.set_xticklabels(main_ax.get_xticklabels(), rotation=45, ha='right')

fig = main_ax.get_figure()
pdf_path = output_base / f"DG_imn_curated_{stamp}.pdf"
png_path = output_base / f"DG_imn_curated_{stamp}.png"
fig.savefig(str(pdf_path), bbox_inches='tight')
fig.savefig(str(png_path), bbox_inches='tight', dpi=150)
print(f"\nSaved: {pdf_path.name}")
print(f"Saved: {png_path.name}")

print(f"\n------ Script completed successfully ------")