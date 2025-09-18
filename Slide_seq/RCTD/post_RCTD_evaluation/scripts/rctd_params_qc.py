import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os

# --- Setup ---
filepath = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/OUT/UPDATED_mouse_gene_orthologs_test/anndata_objects/B03_with_RCTD_mouse.h5ad"
adata = sc.read_h5ad(filepath)
out_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/post_RCTD_evaluation/UPDATED_orthologs_qc"
os.makedirs(out_dir, exist_ok=True)

# --- Compute per-spot and per-gene counts ---
adata.obs["n_UMIs"] = np.array(adata.X.sum(axis=1)).flatten()
adata.var["total_counts"] = np.array(adata.X.sum(axis=0)).flatten()

total_sum = adata.var["total_counts"].sum()
adata.var["frac_total"] = adata.var["total_counts"] / total_sum

# --- 1. UMI distribution ---
plt.figure(figsize=(6,4))
plt.hist(adata.obs["n_UMIs"], bins=100, color="steelblue")
plt.axvline(100, color="red", linestyle="--", label="UMI_min=100")
plt.xlabel("UMIs per spot")
plt.ylabel("Number of spots")
plt.title("Spot UMI distribution")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "UMI_distribution.png"))
plt.close()

# --- 2. Gene fraction distribution ---
plt.figure(figsize=(6,4))
plt.hist(adata.var["frac_total"], bins=200, color="darkgreen")
plt.axvline(0.000125, color="red", linestyle="--", label="gene_cutoff=0.000125")
plt.axvline(2e-4, color="orange", linestyle="--", label="gene_cutoff_reg=2e-04")
plt.xscale("log")
plt.xlabel("Fraction of total counts per gene (log scale)")
plt.ylabel("Number of genes")
plt.title("Gene abundance distribution")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Gene_fraction_distribution.png"))
plt.close()

# --- 3. Cumulative distribution of gene fractions ---
sorted_fracs = np.sort(adata.var["frac_total"])
cumfrac = np.cumsum(sorted_fracs)

plt.figure(figsize=(6,4))
plt.plot(sorted_fracs, cumfrac/cumfrac[-1], color="purple")
plt.axvline(0.000125, color="red", linestyle="--", label="gene_cutoff")
plt.axvline(2e-4, color="orange", linestyle="--", label="gene_cutoff_reg")
plt.xscale("log")
plt.xlabel("Fraction of total counts per gene (log scale)")
plt.ylabel("Cumulative fraction of counts")
plt.title("Cumulative distribution of gene fractions")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Gene_fraction_cumulative.png"))
plt.close()
