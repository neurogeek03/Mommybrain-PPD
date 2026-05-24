import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad(
    "/scratch/mfafouti/Mommybrain-PPD/Slide_tags/Oligo_states/data/"
    "PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
)

sc.pl.umap(
    adata,
    color="treatment",
    palette={"OIL": "#4C9BE8", "CORT": "#E87B4C"},
    title="UMAP — Treatment (OIL vs CORT)",
    frameon=False,
    show=False,
)

out_path = "/scratch/mfafouti/Mommybrain-PPD/Slide_tags/Oligo_states/out/umap_all_cells_treatment.png"
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved to {out_path}")
