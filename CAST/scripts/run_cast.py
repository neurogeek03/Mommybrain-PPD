import os
import torch
import numpy as np
import anndata as ad
import scanpy as sc
import CAST
import warnings
from pathlib import Path
import json
from datetime import datetime
warnings.filterwarnings("ignore")

stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
print(f'------ Script started at {stamp} ------')

print('imports done')

work_dir = Path.cwd().parents[0]
print(f'Current working directroy {work_dir}')
data_path = work_dir / 'data'
query_path = data_path / 'query' / 'rat' 
output_base = work_dir / 'out'
input_path =  output_base / 'cast_mark_input.json' 
adata_dir = output_base / 'objects'
output_path = output_base / 'CAST_Mark' / f'{stamp}'
output_path.mkdir(exist_ok=True, parents=True)

embed_dict_path = output_base / 'CAST_Mark' / '20251125_222722' / 'demo_embed_dict.pt'

# ============================= PROCESS INPUT DICT =============================
print('loading input .json ....')
with open(input_path, "r") as f: input_json = json.load(f)
print('JSON read and is of type:')
print(type(input_json))

# extracting the 2 dicts 
coords_dict = input_json["coords_dict"]

# for key in coords_dict: coords_dict[key] = np.asarray(coords_dict[key])
    
# coords_dict[key] = np.asarray(coords_dict[key]).reshape(-1, 2) 

for key in coords_dict:
    coords_dict[key] = np.asarray(coords_dict[key])
    coords_dict[key] = coords_dict[key].reshape(-1, 2)

exp_dict = input_json["exp_dict"]

samples = ["B01", "B42"]

# print('creating embeddings ....')
# embed_dict = CAST.CAST_MARK(coords_dict,exp_dict,output_path) # coord_dict = coords_raw 

print('loading embed dict ...')
embed_dict = torch.load(embed_dict_path)

print('running k-means ....')
CAST.kmeans_plot_multiple(embed_dict,samples,coords_dict,'demo1',output_path,k=30,dot_size = 10,minibatch=False)