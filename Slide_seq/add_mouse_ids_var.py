import os
import anndata as ad
import pandas as pd

# === Configuration ===
project_folder = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename"
h5ad_dir = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Spatial/anndata_objects"           # change this
ortholog_csv = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/rat_mouse_genes.csv"    # change this
output_dir = os.path.join(project_folder, "with_mouse_orthologs")
os.makedirs(output_dir, exist_ok=True)

# === Load mapping: rat_gene_id -> mouse_gene_stable_ID ===
ortholog_df = pd.read_csv(ortholog_csv)
ortholog_df = ortholog_df.drop_duplicates("rat_ID")  # in case of duplicates
rat_to_mouse = ortholog_df.set_index("rat_ID")["mouse_gene_stable_ID"]

# === Loop through .h5ad files ===
for fname in os.listdir(h5ad_dir):
    if not fname.endswith(".h5ad"):
        continue

    fpath = os.path.join(h5ad_dir, fname)
    print(f"🔄 Processing {fname}...")

    # Load h5ad
    adata = ad.read_h5ad(fpath)

    # Map rat gene IDs (assumed to be adata.var.index) to mouse orthologs
    rat_ids = adata.var.index
    adata.var["mouse_gene_stable_ID"] = rat_ids.map(rat_to_mouse)

    print(adata.var.head())

    # Save updated h5ad
    out_path = os.path.join(output_dir, fname)
    adata.write(out_path)
    print(f"✅ Saved updated file to: {out_path}")
