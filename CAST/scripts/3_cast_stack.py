import CAST
import torch
import numpy as np
import warnings
from pathlib import Path
import json
from datetime import datetime
warnings.filterwarnings("ignore")

stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
print(f'------ Script started at {stamp} ------')

print('imports done')

# paths 
work_dir = Path.cwd().parents[0]
print(f'Current working directroy {work_dir}')

output_base = work_dir / 'out'
input_path =  output_base / 'cast_mark_input.json'
output_path = output_base / 'CAST_Stack' / f'{stamp}'
output_path.mkdir(exist_ok=True, parents=True)

embed_dict_path = output_base / 'CAST_Mark' / '20251125_222722' / 'demo_embed_dict.pt'

# ============================= loading previous data =============================
# coords
print('loading input .json ....')
with open(input_path, "r") as f: input_json = json.load(f)
coords_dict = input_json["coords_dict"]
for key in coords_dict: coords_dict[key] = np.asarray(coords_dict[key]).reshape(-1, 2)
# for key in coords_dict:
#     coords_dict[key] = np.asarray(coords_dict[key])
#     coords_dict[key] = coords_dict[key].reshape(-1, 2)
coords_dict = input_json["coords_dict"]

# coords_dict = torch.load("coords_dict.pt")
coords_to_save = {k: torch.tensor(v) for k, v in coords_dict.items()}
coords_file_path = output_base / "coords_dict.pt"
torch.save(coords_to_save, coords_file_path)

# embeddings
embed_dict = torch.load(embed_dict_path,map_location='cpu')
# graph list
graph_list = ['B01','B42'] # [query_sample, reference_sample]
query_sample = graph_list[0]

# ============================= params =============================
params_dist = CAST.reg_params(dataname = query_sample,
                            gpu = 0 if torch.cuda.is_available() else -1, 
                            diff_step = 5,
                            #### Affine parameters
                            iterations=500,
                            dist_penalty1=0,
                            bleeding=500,
                            d_list = [3,2,1,1/2,1/3],
                            attention_params = [None,3,1,0], 
                            #### FFD parameters    
                            dist_penalty2 = [0],
                            alpha_basis_bs = [500],
                            meshsize = [8],
                            iterations_bs = [400],
                            attention_params_bs = [[None,3,1,0]],
                            mesh_weight = [None])
params_dist.alpha_basis = torch.Tensor([1/1000,1/1000,1/50,5,5]).reshape(5,1).to(params_dist.device)

# ============================= running stack =============================
coords_final = CAST.CAST_STACK(coords_dict,embed_dict,output_path,graph_list,params_dist)