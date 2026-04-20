"""
Fix gene symbols in adata.var["gene_symbol"] by parsing the "name" column.

Format in "name" column: {gene_name}-{ENSEMBL_ID}
  - If gene_name is "." or empty → fall back to Ensembl ID from var_names
  - Otherwise use gene_name

Saves a new h5ad with updated gene_symbol column.
"""

import re
import scanpy as sc

import os

INPUT  = "/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/data/All_RCTD_types_singlet_score_0_slide_seq_15.h5ad"
OUTPUT = "/scratch/mfafouti/Mommybrain/Slide_seq/Spatial/plot_genes/data/All_RCTD_types_singlet_score_0_slide_seq_15_named.h5ad"

os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)


def parse_gene_name(name_val, ensembl_id):
    """Extract gene name from 'name' column entry; fall back to Ensembl ID if absent."""
    match = re.match(r"^(.*)-ENS[A-Z]+\d+$", str(name_val))
    if match:
        gene_name = match.group(1).strip()
        if gene_name and gene_name != ".":
            return gene_name
    return ensembl_id


print("Loading ...")
adata = sc.read_h5ad(INPUT)

print("Parsing gene names from 'name' column ...")
new_symbols = [
    parse_gene_name(name_val, ensembl_id)
    for name_val, ensembl_id in zip(adata.var["name"], adata.var_names)
]

# Summary
n_proper   = sum(1 for s, e in zip(new_symbols, adata.var_names) if s != e)
n_fallback = sum(1 for s, e in zip(new_symbols, adata.var_names) if s == e)
print(f"  Genes with a proper name:        {n_proper}")
print(f"  Genes falling back to Ensembl ID: {n_fallback}")
print(f"  Examples: {list(zip(adata.var_names[:5].tolist(), new_symbols[:5]))}")

adata.var["gene_symbol"] = new_symbols

print(f"Saving to {OUTPUT} ...")
adata.write_h5ad(OUTPUT)
print("Done.")
