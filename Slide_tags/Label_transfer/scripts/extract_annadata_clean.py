import sys
from pathlib import Path

# Add the parent folder (one level up) so Python can find Label_transfer
base_dir = Path(__file__).resolve().parents[3]
print(base_dir)
sys.path.append(str(base_dir / "functions"))

from convert_to_seurat import export_anndata_to_mtx
print('ok!')

# make output dir 
project_dir = Path.cwd().parents[0]
print([project_dir])
current_dir = project_dir / 'query'

export_anndata_to_mtx(
    input_path= current_dir / f'merged_filtered_129084_mincells_10_in_2_samples_slide_tags.h5ad'
)
