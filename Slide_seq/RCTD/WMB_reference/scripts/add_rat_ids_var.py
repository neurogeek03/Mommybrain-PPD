import os
import anndata as ad
import pandas as pd

# === Configuration ===
project_folder = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference_coronal"
h5ad_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/data"          
ortholog_csv = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/rat_mouse_genes.csv"   
output_dir = os.path.join(project_folder, "ABC_regions_mouse_with_rat_genes_h5ads")
os.makedirs(output_dir, exist_ok=True)

# === Load mapping: rat_gene_id -> mouse_gene_stable_ID ===
ortholog_df = pd.read_csv(ortholog_csv)
ortholog_df = ortholog_df.drop_duplicates("mouse_gene_stable_ID")  # in case of duplicates
mouse_to_rat = ortholog_df.set_index("mouse_gene_stable_ID")["rat_ID"]

# === Loop through .h5ad files ===
for fname in os.listdir(h5ad_dir):
    if not fname.endswith(".h5ad"):
        continue

    fpath = os.path.join(h5ad_dir, fname)
    print(f"🔄 Processing {fname}...")

    # Load h5ad
    adata = ad.read_h5ad(fpath)

    # Map rat gene IDs (assumed to be adata.var.index) to mouse orthologs
    mouse_ids = adata.var.index # column: "gene_identifier"
    adata.var["rat_ID"] = mouse_ids.map(mouse_to_rat)

    # Count how many mapped
    total_genes = len(mouse_ids)
    mapped_genes = adata.var["rat_ID"].notna().sum()
    print(f"   🧬 Mapped {mapped_genes:,} / {total_genes:,} genes ({mapped_genes/total_genes:.1%})")

    print(adata.var.sample(10, random_state=42))

    # Save updated h5ad
    out_path = os.path.join(output_dir, fname)
    adata.write(out_path)
    print(f"✅ Saved updated file to: {out_path}")
