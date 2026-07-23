import scanpy as sc
import matplotlib.pyplot as plt

DATA = "data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
OUT = "out"
OLIGO_SUBCLASS = "327 Oligo NN"
RESOLUTION = 0.0
N_HVG = 3000
N_PCS = 4
N_NEIGHBORS = 15

adata = sc.read_h5ad(DATA)
oligo = adata[adata.obs["subclass_name"] == OLIGO_SUBCLASS].copy()
print(f"Oligo NN cells: {oligo.n_obs}")
print(oligo.obs["treatment"].value_counts())

sc.settings.figdir = OUT

# .X is already scaled — go straight to PCA
sc.tl.pca(oligo, n_comps=N_PCS)
sc.pl.pca_variance_ratio(oligo, n_pcs=N_PCS, log=True,
                          save=f"_oligo_elbow_res{RESOLUTION}.png", show=False)
sc.pp.neighbors(oligo, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
sc.tl.umap(oligo)
sc.tl.leiden(oligo, resolution=RESOLUTION, key_added="leiden")

# UMAP: treatment
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.umap(oligo, color="treatment", palette={"OIL": "#4C72B0", "CORT": "#DD8452"},
           title="Treatment", ax=axes[0], show=False)
sc.pl.umap(oligo, color="leiden", title=f"Leiden (res={RESOLUTION})",
           ax=axes[1], show=False, legend_loc="on data")
plt.tight_layout()
fig.savefig(f"{OUT}/oligo_umap_treatment_leiden_res{RESOLUTION}.png", dpi=150, bbox_inches="tight")
plt.close()

# UMAP: treatment split side by side
sc.pl.umap(oligo, color="treatment", groups=["OIL", "CORT"],
           palette={"OIL": "#4C72B0", "CORT": "#DD8452"},
           title="Treatment", save=f"_oligo_treatment_res{RESOLUTION}.png", show=False)

oligo.write_h5ad(f"{OUT}/oligo_NN_reclustered_res{RESOLUTION}.h5ad")
print("Done. Outputs written to out/")
