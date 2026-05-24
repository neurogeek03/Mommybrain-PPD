import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path


def get_broad_class(subclass_name):
    if subclass_name.endswith('Glut'):
        return 'Glutamatergic'
    elif subclass_name.endswith('Gaba'):
        return 'GABAergic'
    elif subclass_name.endswith('NN'):
        return 'Non-neuronal'
    return 'Unknown'


def load_cell_types_from_edger(folder):
    """
    Parse subclass names from EdgeR results filenames.
    Expected pattern: NNN_SubclassName_edgeR_results.tsv
    Returns list of subclass name strings.
    """
    folder = Path(folder)
    names = []
    for f in sorted(folder.glob('*_edgeR_results.tsv')):
        stem = f.stem.replace('_edgeR_results', '')
        names.append(stem)
    return names


def check_log_normalized(adata):
    signal1 = 'log1p' in adata.uns
    if sp.issparse(adata.X):
        max_val = adata.X.max()
    else:
        max_val = float(np.max(adata.X)) if adata.X.size > 0 else 0.0
    signal2 = max_val < 15
    print("--- Log-Normalization Diagnostic ---")
    print(f"  'log1p' in adata.uns : {signal1}")
    print(f"  max(adata.X) < 15   : {signal2}  (max={max_val:.2f})")
    print(f"  adata.raw is not None: {adata.raw is not None}")
    print(f"  Likely log-normalized: {signal1 or signal2}\n")
    return signal1 or signal2


def compute_dotplot_data(adata, marker_genes, group_col, subclass_order):
    """
    Compute per-(gene, subclass) mean expression, pct expressing, and specificity.

    Returns a DataFrame with columns:
        subclass, gene, mean_expr, pct_expr, specificity, broad_class, subclass_order
    Genes absent from adata.var_names are skipped with a warning.
    """
    present = [g for g in marker_genes if g in adata.var_names]
    missing = [g for g in marker_genes if g not in adata.var_names]
    if missing:
        print(f"WARNING: {len(missing)} markers not in data (skipped): {missing}")

    adata_m = adata[:, present]
    groups = list(subclass_order)

    pct_df = pd.DataFrame(index=present, columns=groups, dtype=float)
    mean_df = pd.DataFrame(index=present, columns=groups, dtype=float)

    for group in groups:
        mask = adata_m.obs[group_col] == group
        n = int(mask.sum())
        if n == 0:
            pct_df[group] = 0.0
            mean_df[group] = 0.0
            continue
        X = adata_m[mask].X
        arr = X.toarray() if sp.issparse(X) else np.asarray(X)
        pct_df[group] = (arr > 0).mean(axis=0)
        mean_df[group] = arr.mean(axis=0)

    spec_df = pd.DataFrame(index=present, columns=groups, dtype=float)
    for group in groups:
        other_max = pct_df.drop(columns=group).max(axis=1)
        spec_df[group] = pct_df[group] - other_max

    rows = []
    for i, group in enumerate(subclass_order):
        for gene in present:
            rows.append({
                'subclass': group,
                'gene': gene,
                'mean_expr': float(mean_df.at[gene, group]),
                'pct_expr': float(pct_df.at[gene, group]),
                'specificity': float(spec_df.at[gene, group]),
                'broad_class': get_broad_class(group),
                'subclass_order': i,
            })

    return pd.DataFrame(rows)
