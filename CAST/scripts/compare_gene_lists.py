import scanpy as sc
from pathlib import Path
import pandas as pd

project_path = Path.cwd()
ref_dir = project_path / 'data' / 'reference' / 'Xiaowei' / 'ABCA_1'
print(ref_dir)
own_dir = project_path / 'data' /'query'/'mouse'

ref_object = ref_dir / 'Zhuang-ABCA-1-raw.h5ad'
query_object = own_dir / 'B01_collapsed_mouse_genes.h5ad'

# read in data
ref_ad = sc.read_h5ad(ref_object)
query_ad = sc.read_h5ad(query_object)

genes_ref = set(ref_ad.var_names)
genes_query = set(query_ad.var_names)

shared_genes = genes_ref.intersection(genes_query)

print(f"Total genes in reference: {len(genes_ref)}")
print(f"Total genes in query: {len(genes_query)}")
print(f"Shared genes: {len(shared_genes)}")

print('test2!')