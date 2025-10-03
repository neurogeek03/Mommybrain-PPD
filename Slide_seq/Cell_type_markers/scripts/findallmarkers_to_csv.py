import scanpy as sc
import pandas as pd

# Load your AnnData object
adata = sc.read_h5ad("/scratch/mfafouti/Mommybrain/Slide_seq/Cell_type_markers/data/filtered_10_subclass_sample_umap_cleaned_mouse_RCTD_slideseq_singlets_15samples.h5ad")

# Make sure cluster IDs are in adata.obs["clusters"]
# (replace "clusters" with your actual clustering column)
sc.tl.rank_genes_groups(
    adata,
    groupby="clusters",
    method="wilcoxon",   # alternatives: 't-test', 'logreg'
    n_genes=adata.shape[1]
)

# Extract results into a dataframe
markers = pd.DataFrame()
groups = adata.obs["clusters"].unique().tolist()

for group in groups:
    result = sc.get.rank_genes_groups_df(adata, group=group)
    result["cluster"] = group
    markers = pd.concat([markers, result], ignore_index=True)

# Optional: add a "pct_diff" like in Seurat
# proportion of cells in cluster expressing gene minus other clusters
sc.tl.dendrogram(adata, groupby="clusters")  # ensures categories ordered
for group in groups:
    mask = adata.obs["clusters"] == group
    pct_in = (adata[mask].X > 0).sum(axis=0).A1 / mask.sum()
    pct_out = (adata[~mask].X > 0).sum(axis=0).A1 / (~mask).sum()
    markers.loc[markers["cluster"] == group, "pct.1"] = pct_in[markers.loc[markers["cluster"] == group, "names"].map(lambda g: adata.var_names.get_loc(g))]
    markers.loc[markers["cluster"] == group, "pct.2"] = pct_out[markers.loc[markers["cluster"] == group, "names"].map(lambda g: adata.var_names.get_loc(g))]
    markers.loc[markers["cluster"] == group, "pct_diff"] = markers.loc[markers["cluster"] == group, "pct.1"] - markers.loc[markers["cluster"] == group, "pct.2"]

# Save to CSV
markers.to_csv("all_markers.csv", index=False)
