# Marker Gene Selection and Dotplot CSV Export

How `create_marker_csv.py` selects marker genes and builds the output CSV
consumed by `plot_markers.py`.

## Pipeline overview

```
h5ad  -->  preprocess  -->  rank genes  -->  specificity scoring  -->  select top markers  -->  CSV export
```

---

## 1. Preprocessing

- **Normalize**: `sc.pp.normalize_total` (library-size normalization to median counts).
- **Log-transform**: `sc.pp.log1p`.
- Gene symbols from `adata.var['name']` replace Ensembl IDs in `var_names`.

All downstream computations (ranking, mean expression, percent expressing) operate
on this log-normalized matrix.

## 2. Subclass filtering

- If `subclasses/shared_subclasses.csv` exists, only those subclasses are kept.
  Otherwise all subclasses in `adata.obs['RCTD_first_type_rat']` are used.
- Subclasses are sorted for display: first by broad class
  (Glutamatergic > GABAergic > Non-neuronal > Unknown), then alphabetically within
  each broad class. Broad class is inferred from the subclass name suffix
  (`Glut`, `Gaba`, `NN`).

## 3. Differential expression (candidate pool)

- `sc.tl.rank_genes_groups` with Wilcoxon rank-sum test, one-vs-rest per subclass,
  on the log-normalized matrix (`use_raw=False`).
- For each subclass, the top `RANK_POOL_SIZE` (default 50) genes by adjusted p-value
  are collected as the candidate pool.

## 4. Specificity scoring

For every candidate gene *g* and subclass *s*:

1. **Percent expressing** (`pct_expr`): fraction of cells in subclass *s* where
   expression of *g* > 0.
2. **Specificity score**: `pct_expr(g, s) - max over all other subclasses(pct_expr(g, s'))`.

A positive specificity means the gene is expressed in a larger fraction of cells in
that subclass than in any other single subclass. Higher is more specific.

The full specificity matrix (all candidates x all subclasses) is saved as
`out/candidate_specificity_scores_*.csv`.

## 5. Top marker selection

For each subclass (in display order):

1. Take its candidate pool from step 3.
2. Sort those candidates by their specificity score for that subclass (descending).
3. Optionally drop genes below `MIN_SPEC_SCORE` (default 0, i.e. no filter).
4. Keep the top `TOP_N_GENES` (default 2).

Genes are assembled into a flat list following subclass display order, producing a
diagonal pattern on the dotplot: each subclass's markers appear in sequence.

If `CURATED_MARKERS` is populated instead, those genes are used directly and the
automatic selection above is skipped.

## 6. Dotplot data CSV export

The script subsets `adata` to the final marker gene list and computes, for every
(gene, subclass) pair:

| Column | Description |
|--------|-------------|
| `subclass` | Subclass label |
| `gene` | Gene symbol |
| `mean_expr` | Mean log-normalized expression across cells in that subclass |
| `pct_expr` | Fraction of cells in that subclass expressing the gene (> 0) |
| `specificity` | Specificity score from step 4 |
| `broad_class` | Glutamatergic / GABAergic / Non-neuronal / Unknown |
| `subclass_order` | Integer rank preserving the display order |

Output: `out/marker_dotplot_data_<timestamp>.csv`.

This CSV is the sole input to `plot_markers.py`, which reads it and renders the
dotplot without needing scanpy or the h5ad file.

## Key parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| `RANK_POOL_SIZE` | 50 | Candidate genes per subclass from Wilcoxon ranking |
| `TOP_N_GENES` | 2 | Final markers per subclass after specificity re-ranking |
| `MIN_SPEC_SCORE` | 0.0 | Floor on specificity (0 = keep all) |
| `CURATED_MARKERS` | `{}` | If non-empty, bypasses automatic selection entirely |
