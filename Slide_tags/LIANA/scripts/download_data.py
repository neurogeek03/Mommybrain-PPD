import liana as li
import decoupler as dc
from pathlib import Path
import pandas as pd
import scanpy as sc
import glob

net = dc.op.collectri(organism='human', remove_complexes=False, license='academic', verbose=False)