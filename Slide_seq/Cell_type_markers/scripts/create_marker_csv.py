import os
import scanpy as sc
import pandas as pd
import numpy as np

out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/Cell_type_markers/out"
ad_file = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/anndata_objects/delta_5_umi30_subclass_B03_with_RCTD_mouse.h5ad"

# --- Read data ---
adata = sc.read_h5ad(ad_file)
print(adata)

#subsetting to singlets
adata = adata[adata.obs["RCTD_spot_class_mouse"] == "singlet"].copy()
print(f"Remaining cells: {adata.n_obs}, genes: {adata.n_vars}")

adata.var_names = adata.var["name"]
adata.var_names_make_unique()

group_col = "RCTD_first_type_mouse"
groups = adata.obs[group_col].unique()

# --- Specificity calculation ---
binary_expr = adata.raw.X > 0 if adata.raw else adata.X > 0
pct_expr = pd.DataFrame(index=adata.var_names, columns=groups)

for group in groups:
    idx = adata.obs[group_col] == group
    group_expr = binary_expr[idx.values, :]
    pct_expr[group] = np.asarray(group_expr.sum(axis=0)).flatten() / idx.sum()

specificity = pd.DataFrame(index=adata.var_names, columns=groups)
for group in groups:
    specificity[group] = pct_expr[group] - pct_expr.drop(columns=group).max(axis=1)

best_class = specificity.idxmax(axis=1)
best_scores = specificity.max(axis=1)
gene_best = pd.DataFrame({
    'gene': specificity.index,
    'best_subclass': best_class,
    'specificity': best_scores
})

# Save specificity tables
gene_best.to_csv(os.path.join(out_dir, "gene_best_subclass_specificity.csv"))
specificity.to_csv(os.path.join(out_dir, "subclass_gene_specificity.csv"))

# --- Pick top markers for plotting ---
top_genes = {}
for group in groups:
    top_genes[group] = specificity[group].sort_values(ascending=False).head(2).index.tolist()
genes_for_plot = list({gene for genes in top_genes.values() for gene in genes})

# --- Dotplot-style summary ---
results = []
for group in groups:
    cells_in_group = adata.obs[group_col] == group
    if cells_in_group.sum() == 0:  # skip empty groups
        continue

    subset = adata[cells_in_group]

    for gene in genes_for_plot:
        if gene not in adata.var_names:
            continue

        expr = subset[:, gene].X
        if not isinstance(expr, np.ndarray):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()

        if expr.size == 0:  # skip empty slices
            continue

        results.append({
            "group": group,
            "gene": gene,
            "mean_expression": float(expr.mean()),
            "pct_expressing": float((expr > 0).mean())
        })


df = pd.DataFrame(results)
df.to_csv(os.path.join(out_dir, "dotplot_data.csv"), index=False)

