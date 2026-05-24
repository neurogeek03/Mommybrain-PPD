import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

adata = sc.read_h5ad("data/All_RCTD_types_B03_B14_filtered.h5ad")

# rename var to symbols
adata.var_names = adata.var["gene_symbol"].values
adata.var_names_make_unique()

# ── subset: Astroependymal singlets only ──────────────────────────────────────
mask = (
    (adata.obs["RCTD_first_type_rat"] == "321_Astroependymal_NN") &
    (adata.obs["RCTD_spot_class_rat"] == "singlet")
)
sub = adata[mask].copy()
print(f"Astroependymal singlets: {sub.n_obs}")

# ── preprocessing ─────────────────────────────────────────────────────────────
sc.pp.normalize_total(sub, target_sum=1e4)
sc.pp.log1p(sub)
sc.pp.highly_variable_genes(sub, n_top_genes=3000, flavor="seurat")
sub = sub[:, sub.var["highly_variable"]].copy()
sc.pp.scale(sub, max_value=10)
sc.tl.pca(sub, n_comps=30)

# ── neighbors + UMAP + Leiden ─────────────────────────────────────────────────
sc.pp.neighbors(sub, n_pcs=20)
sc.tl.umap(sub)
sc.tl.leiden(sub, resolution=0.3, key_added="leiden", flavor="igraph", directed=False, n_iterations=2)

print("Leiden cluster sizes:")
print(sub.obs["leiden"].value_counts().sort_index())

# ── UMAP plots ────────────────────────────────────────────────────────────────
sc.settings.figdir = "plots/"

fig, axes = plt.subplots(1, 1, figsize=(6, 5))
sc.pl.umap(sub, color="leiden", ax=axes, show=False, title="Leiden clusters (res=0.3)")
plt.tight_layout()
plt.savefig("plots/astroepen_umap_leiden.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved plots/astroepen_umap_leiden.png")

# ── marker gene UMAP panels ───────────────────────────────────────────────────
all_markers = ["Aqp4","Cd38","Agt","Itih3","Rfx4","Dbx2","Prdm16","Nr2f1",
               "Gja1","Lhx2","Gpc5","Sox9","Meis2",
               "Gfap","Slit2","Slc7a11","Zic4","Zfhx4","Prrx1"]
present = [g for g in all_markers if g in sub.var_names]
print("Markers present for UMAP:", present)

n = len(present)
ncols = 5
nrows = -(-n // ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3))
axes = axes.flatten()
for i, g in enumerate(present):
    sc.pl.umap(sub, color=g, ax=axes[i], show=False, title=g, colorbar_loc=None)
for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)
plt.suptitle("Marker gene expression — Astroependymal reclustered", y=1.01, fontsize=12)
plt.tight_layout()
plt.savefig("plots/astroepen_umap_markers.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved plots/astroepen_umap_markers.png")

# ── dotplot by new cluster ────────────────────────────────────────────────────
# re-read raw counts for dotplot (scale ruins expression values)
adata2 = sc.read_h5ad("data/All_RCTD_types_B03_B14_filtered.h5ad")
adata2.var_names = adata2.var["gene_symbol"].values
adata2.var_names_make_unique()
sub2 = adata2[sub.obs_names].copy()
sub2.obs["leiden"] = sub.obs["leiden"]
sc.pp.normalize_total(sub2, target_sum=1e4)
sc.pp.log1p(sub2)

markers_plot = {
    "Astro_NT — subclass": ["Aqp4","Cd38","Agt","Itih3"],
    "Astro_NT — TF":       ["Rfx4","Dbx2","Prdm16","Nr2f1"],
    "Astro_TE — subclass": ["Gja1","Lhx2","Gpc5","Nr2f1"],
    "Astro_TE — TF":       ["Rfx4","Lhx2","Sox9","Nr2f1","Meis2"],
    "Astroepen — subclass":["Gfap","Slit2","Slc7a11","Zic4"],
    "Astroepen — TF":      ["Rfx4","Zic4","Zfhx4","Nr2f1","Prrx1"],
}
markers_plot = {k: [g for g in v if g in sub2.var_names] for k, v in markers_plot.items()}

dp = sc.pl.dotplot(
    sub2,
    var_names=markers_plot,
    groupby="leiden",
    standard_scale="var",
    dot_max=0.8,
    figsize=(18, 5),
    title="Marker genes by Leiden cluster — Astroependymal",
    show=False,
    return_fig=True,
)
dp.savefig("plots/astroepen_cluster_dotplot.png", dpi=150, bbox_inches="tight")
print("Saved plots/astroepen_cluster_dotplot.png")
