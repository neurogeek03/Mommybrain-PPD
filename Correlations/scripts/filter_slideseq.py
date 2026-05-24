import scanpy as sc

INPUT  = "data/All_RCTD_types_singlet_score_0_slide_seq_15_named.h5ad"
OUTPUT = "data/All_RCTD_types_B03_B14_filtered.h5ad"

adata = sc.read_h5ad(INPUT)
print(f"Loaded: {adata.shape}")

# Sample filter
adata = adata[adata.obs["sample"].isin(["B03", "B14"])].copy()
print(f"After sample filter (B03, B14): {adata.shape}")

is_nn = adata.obs["RCTD_first_type_rat"].str.endswith("NN", na=False)

# Non-neurons: singlet class + score > 330
mask_nn = is_nn & (adata.obs["RCTD_spot_class_rat"] == "singlet") & (adata.obs["RCTD_singlet_score_rat"] > 330)

# Neurons / other: score > 330 only
mask_other = ~is_nn & (adata.obs["RCTD_singlet_score_rat"] > 330)

adata = adata[mask_nn | mask_other].copy()
print(f"After quality filter: {adata.shape}")
print(f"  NN spots:    {mask_nn.sum()}")
print(f"  Other spots: {mask_other.sum()}")

# Merge Astro-NT, Astro-TE, and Astroependymal into a single "Astrocytes" label
# Color: #AD5CCC (Astro-NT NN from cluster_annotation_term.csv)
ASTRO_TYPES = {"318_Astro_NT_NN", "319_Astro_TE_NN", "321_Astroependymal_NN"}
ASTRO_COLOR = "#AD5CCC"

adata.obs["cell_type"] = adata.obs["RCTD_first_type_rat"].astype(str).copy()
is_astro = adata.obs["cell_type"].isin(ASTRO_TYPES)
adata.obs.loc[is_astro, "cell_type"] = "Astrocytes"
print(f"Merged {is_astro.sum()} spots into 'Astrocytes'")

# Log-normalize to match slide-tags scale (log(CP10k + 1))
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("Log-normalized (log(CP10k+1))")

# Store color mapping in uns for downstream use
color_map = {"Astrocytes": ASTRO_COLOR}
adata.uns["cell_type_colors"] = color_map

adata.write_h5ad(OUTPUT)
print(f"Saved to {OUTPUT}")