import scanpy as sc

h5ad_path = "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/classifier/out/78689_label_transfer_scanorama.h5ad"
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/classifier/out"
fig_name = "78689_umap_predicted_label.png"

adata = sc.read_h5ad(h5ad_path)

if "X_umap" not in adata.obsm:
    sc.pp.neighbors(adata, use_rep="X_scanorama")
    sc.tl.umap(adata)

sc.settings.figdir = out_dir

sc.pl.umap(
    adata,
    color="predicted_label",
    legend_loc="on data",
    frameon=False,
    save=f"_{fig_name}",
    show=False,
)
