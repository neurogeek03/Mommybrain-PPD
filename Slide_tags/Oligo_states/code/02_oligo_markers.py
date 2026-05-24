import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

OLIGO_H5AD = "out/oligo_NN_reclustered_res0.1.h5ad"
RAW_H5AD = "data/DE_after_mt_filter_108123_mincells_10_in_2_samples_slide_tags.h5ad"
OUT = "out"
CLUSTERS_OF_INTEREST = ["0", "4"]
OLIGO_SUBCLASS = "327 Oligo NN"

sc.settings.figdir = OUT

# Load reclustered object to get leiden labels
oligo_scaled = sc.read_h5ad(OLIGO_H5AD)

# Load raw counts, subset to oligo cells present in the reclustered object
raw = sc.read_h5ad(RAW_H5AD)
raw.var_names = raw.var["gene_symbol"].astype(str)
raw.var_names_make_unique()
oligo = raw[raw.obs["subclass_name"] == OLIGO_SUBCLASS].copy()

# Transfer leiden labels (align on cell barcodes)
oligo.obs["leiden"] = oligo_scaled.obs["leiden"].reindex(oligo.obs_names)
oligo = oligo[oligo.obs["leiden"].notna()].copy()
print(oligo.obs["leiden"].value_counts())

# Normalize and log1p for proper DE
sc.pp.normalize_total(oligo, target_sum=1e4)
sc.pp.log1p(oligo)

# Marker genes: each cluster vs all others
sc.tl.rank_genes_groups(oligo, groupby="leiden", method="wilcoxon", key_added="markers")

# Export full marker table
dfs = []
for cluster in oligo.obs["leiden"].unique():
    df = sc.get.rank_genes_groups_df(oligo, group=cluster, key="markers")
    df.insert(0, "cluster", cluster)
    dfs.append(df)
markers_df = pd.concat(dfs).reset_index(drop=True)
markers_df.to_csv(f"{OUT}/oligo_markers_all_clusters.csv", index=False)

# Focused table: top 50 markers for clusters of interest
for c in CLUSTERS_OF_INTEREST:
    df = sc.get.rank_genes_groups_df(oligo, group=c, key="markers").head(50)
    df.to_csv(f"{OUT}/oligo_markers_cluster{c}.csv", index=False)
    print(f"\nTop 20 markers for cluster {c}:")
    print(df[["names", "scores", "logfoldchanges", "pvals_adj"]].head(20).to_string(index=False))

# Dotplot: top 5 markers per cluster
sc.pl.rank_genes_groups_dotplot(
    oligo, groupby="leiden", key="markers", n_genes=5,
    min_logfoldchange=None,
    save="_oligo_markers_dotplot.png", show=False
)

# Heatmap of top 10 markers for clusters 0 and 4
top_genes = []
for c in CLUSTERS_OF_INTEREST:
    df = sc.get.rank_genes_groups_df(oligo, group=c, key="markers")
    top_genes += df["names"].head(10).tolist()
top_genes = list(dict.fromkeys(top_genes))  # deduplicate, preserve order

sc.pl.dotplot(
    oligo, var_names=top_genes, groupby="leiden",
    title=f"Top markers for clusters {CLUSTERS_OF_INTEREST}",
    save=f"_oligo_cluster{'_'.join(CLUSTERS_OF_INTEREST)}_markers.png", show=False
)

print("\nDone. Outputs written to out/")
