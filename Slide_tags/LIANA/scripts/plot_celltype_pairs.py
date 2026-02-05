import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import glob

import plotnine as p9
import liana as li

project_path = Path.cwd().parents[0]
lr_res_path = project_path / 'out'/ 'first_try' / 'lr_liana_res.csv'
output_dir = project_path / 'out' / 'celltype_interact'
output_dir.mkdir(exist_ok=True, parents=True)

lr_res = pd.read_csv(lr_res_path)

pd.set_option('display.max_columns', None, index_col=0)
print(lr_res.head())