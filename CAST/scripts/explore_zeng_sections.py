from functions import plot_brains
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path 

project_folder = Path.cwd() # always run from project folder 
metadata_file = project_folder / 'data' / 'reference' / 'meta' / 'cell_metadata.csv'
group_column = 'brain_section_label'
color_column = 'cluster_alias'

plot_brains(metadata_file, group_column, 1)
# Example usage:
# plot_groups("my_file.csv", group_col="slice_id")
