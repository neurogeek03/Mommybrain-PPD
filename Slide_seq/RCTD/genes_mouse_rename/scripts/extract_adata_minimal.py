import anndata as ad
import scipy.io
import pandas as pd
import os

project_path = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename"
adata_dir = os.path.join(project_path, "UPDATED_collapsed_adata_with_mouse_orthologs")
out_root = os.path.join(project_path, "UPDATED_rds_files_mouse")
os.makedirs(out_root, exist_ok=True)

def export_minimal_for_r(adata, out_dir, basename):
    os.makedirs(out_dir, exist_ok=True)

    # === Preview first few entries ===
    print(f"\n📦 Preview of exported minimal data for {basename}")
    print("🧬 Genes (first 5):")
    print(adata.var.index[:5].tolist())
    
    print("\n🧫 Barcodes / Cell IDs (first 5):")
    print(adata.obs.index[:5].tolist())
    
    # Save expression matrix
    scipy.io.mmwrite(os.path.join(out_dir, f"{basename}_matrix.mtx"), adata.X)

    # Save genes and cells
    adata.var.index.to_series().to_csv(os.path.join(out_dir, f"{basename}_genes.tsv"), sep="\t", index=False, header=False)
    adata.obs.index.to_series().to_csv(os.path.join(out_dir, f"{basename}_barcodes.tsv"), sep="\t", index=False, header=False)

    print(f"✅ Exported minimal data for {basename}")


# === Loop through all .h5ad files ===
for fname in os.listdir(adata_dir):
    if not fname.endswith(".h5ad"):
        continue

    print(f"🔄 Processing {fname}...")

    fpath = os.path.join(adata_dir, fname)
    adata = ad.read_h5ad(fpath)

    # Remove ".h5ad" and use the rest as the basename
    basename = fname.split("_")[0]

    # Make a subdirectory for each output (optional, or just use out_root)
    out_dir = os.path.join(out_root, basename)
    
    export_minimal_for_r(adata, out_dir=out_dir, basename=basename)