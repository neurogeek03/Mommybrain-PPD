import os
import anndata as ad
import pandas as pd
from pathlib import Path
from datetime import datetime

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# =================== PARAMS ===================
sample ='MT'

# === Configuration ===
project_path = Path.cwd().parents[0]
# subfolder = project_path / 'slideseq_test'
subfolder = Path('/scratch/mfafouti/Mommybrain/Slide_tags/Neuronchat')
output_base = subfolder / 'out'
h5ad_dir = subfolder / 'data' / 'raw'
ortholog_csv = project_path / 'rat_to_mouse_filtered.csv'
output_dir = output_base / 'add_mouse_orthologs'
output_dir.mkdir(exist_ok=True, parents=True)

# === Load mapping: rat_gene_id -> mouse_gene_stable_ID ===
ortholog_df = pd.read_csv(ortholog_csv)
# ortholog_df = ortholog_df.drop_duplicates("rat_ID")  # in case of duplicates
rat_to_mouse = ortholog_df.set_index("Gene stable ID")["Mouse gene stable ID"]

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
    adata.var["mouse_ID"] = rat_ids.map(rat_to_mouse)

    # Count how many mapped
    total_genes = len(rat_ids)
    mapped_genes = adata.var["mouse_ID"].notna().sum()
    print(f"   🧬 Mapped {mapped_genes:,} / {total_genes:,} genes ({mapped_genes/total_genes:.1%})")

    print(adata.var.head())

    # Save updated h5ad
    out_path = os.path.join(output_dir, fname)
    adata.write(out_path)
    print(f"✅ Saved updated file to: {out_path}")
