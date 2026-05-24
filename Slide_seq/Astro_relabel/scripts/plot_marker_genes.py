import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

adata = sc.read_h5ad("data/All_RCTD_types_B03_B14_filtered.h5ad")

# ── build symbol → ensembl ID map ─────────────────────────────────────────────
sym2id = (
    adata.var["gene_symbol"]
    .reset_index()
    .set_index("gene_symbol")["gene_id"]
    .to_dict()
)

# ── marker gene sets (by symbol) ──────────────────────────────────────────────
markers_sym = {
    "Astro_NT — subclass":       ["Aqp4", "Cd38", "Agt", "Itih3"],
    "Astro_NT — TF":             ["Rfx4", "Dbx2", "Prdm16", "Nr2f1"],
    "Astro_TE — subclass":       ["Gja1", "Lhx2", "Gpc5", "Nr2f1"],
    "Astro_TE — TF":             ["Rfx4", "Lhx2", "Sox9", "Nr2f1", "Meis2"],
    "Astroependymal — subclass": ["Gfap", "Slit2", "Slc7a11", "Zic4"],
    "Astroependymal — TF":       ["Rfx4", "Zic4", "Zfhx4", "Nr2f1", "Prrx1"],
}

# convert to Ensembl IDs, keep symbols as display labels
markers_id = {}
for group, syms in markers_sym.items():
    ids = [sym2id[s] for s in syms if s in sym2id]
    markers_id[group] = ids

# rename var_names to gene symbols for readable axis labels
adata.var_names = adata.var["gene_symbol"].values
adata.var_names_make_unique()

# ── subset: three astro types, rat singlets ────────────────────────────────────
target_types = [
    "318_Astro_NT_NN",
    "319_Astro_TE_NN",
    "321_Astroependymal_NN",
]
mask = (
    adata.obs["RCTD_first_type_rat"].isin(target_types) &
    (adata.obs["RCTD_spot_class_rat"] == "singlet")
)
sub = adata[mask].copy()

# clean label
sub.obs["cell_type"] = (
    sub.obs["RCTD_first_type_rat"]
    .str.replace(r"^\d+_", "", regex=True)
    .str.replace("_NN", "", regex=False)
    .astype("category")
)

print("Cell counts:")
print(sub.obs["cell_type"].value_counts())

# normalise
sc.pp.normalize_total(sub, target_sum=1e4)
sc.pp.log1p(sub)

# convert ID-keyed groups to symbol-keyed (now var_names are symbols)
markers_plot = {}
for group, ids in markers_id.items():
    syms_present = [adata.var.loc[i, "gene_symbol"] if i in adata.var.index else None for i in ids]
    # after renaming, var_names are symbols; look up directly
    syms_in_sub = [s for s in markers_sym[group] if s in sub.var_names]
    if syms_in_sub:
        markers_plot[group] = syms_in_sub

print("\nMarker groups for plot:", {k: v for k, v in markers_plot.items()})

# ── dotplot ───────────────────────────────────────────────────────────────────
dp = sc.pl.dotplot(
    sub,
    var_names=markers_plot,
    groupby="cell_type",
    standard_scale="var",
    dot_max=0.8,
    figsize=(18, 4),
    title="Astrocyte subtype marker genes",
    show=False,
    return_fig=True,
)
dp.savefig("plots/astro_marker_dotplot.png", dpi=150, bbox_inches="tight")
print("Saved plots/astro_marker_dotplot.png")

# ── heatmap ───────────────────────────────────────────────────────────────────
all_syms = []
seen = set()
for syms in markers_plot.values():
    for s in syms:
        if s not in seen:
            all_syms.append(s)
            seen.add(s)

sc.pl.heatmap(
    sub,
    var_names=all_syms,
    groupby="cell_type",
    standard_scale="var",
    figsize=(16, 3),
    show=False,
    swap_axes=False,
)
plt.savefig("plots/astro_marker_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved plots/astro_marker_heatmap.png")
