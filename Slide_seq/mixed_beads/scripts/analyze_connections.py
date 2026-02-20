
import pandas as pd
import os
import plotly.graph_objects as go
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2



project_path = Path.cwd().parents[0]
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)
csv_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/OUT/FINAL_RCTD_newgenelist/merged_metadata_csvs"
metadata_file = "/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR/slide_seq_metadata.csv"
output_file = output_base / "merged_slide_seq_metadata.csv"

# Read and process data
df_list = [pd.read_csv(os.path.join(csv_dir, f)).assign(sample=f.split('_')[0]) for f in os.listdir(csv_dir) if f.endswith(".csv")]
merged_df = pd.concat(df_list, ignore_index=True)
metadata_df = pd.read_csv(metadata_file)
final_df = pd.merge(merged_df, metadata_df,     
                    on='sample', how='left').rename(columns={'treatment': 'condition'})
final_df.to_csv(output_file, index=False)

print(f"Successfully merged {len(df_list)} files and saved to {output_file}")
print("Value counts of 'condition' column:", final_df['condition'].value_counts())

# --- Data Splitting and Plotting ---
final_df['condition'] = final_df['condition'].astype(str)
final_df = final_df[final_df['spot_class'] != 'reject'].copy()
final_df = final_df[final_df['spot_class'] != 'doublet_certain'].copy()
final_df = final_df[final_df['spot_class'] != 'doublet_uncertain'].copy()
cort_mask = final_df['condition'].str.contains('CORT', na=False)
df_cort = final_df[cort_mask]
df_other = final_df[~cort_mask]


# --- Venn Diagram for 'first_type' ---
print("Generating Venn Diagram for 'first_type' cell types ---")

# Get unique 'first_type' values from each group
cort_first_types = set(df_cort['first_type'].unique())
other_first_types = set(df_other['first_type'].unique())

# Create the Venn diagram
try:
    plt.figure(figsize=(8, 8))
    venn2([cort_first_types, other_first_types], ('CORT', 'OTHER'))
    venn_file = output_base / "first_type_venn_diagram.png"
    plt.title("Overlap of 'first_type' Cell Types between CORT and OTHER groups")
    plt.savefig(venn_file)
    plt.close() # Close the figure to free memory
    print(f"Successfully created Venn diagram and saved to {venn_file}")
except ImportError:
    print("Could not create Venn diagram ---")
    print("This feature requires the 'matplotlib_venn' package.")
    print("Please install it by running: pip install matplotlib-venn")
    print("------------------------------------")
except Exception as e:
    print(f"An error occurred while creating the Venn diagram ---{e}")


# List unique and common cell types
unique_to_cort = cort_first_types - other_first_types
unique_to_other = other_first_types - cort_first_types
common_types = cort_first_types.intersection(other_first_types)

print("--- 'first_type' Cell Type Analysis ---")
print(f"Total unique 'first_type' in CORT: {len(cort_first_types)}")
print(f"Total unique 'first_type' in OTHER: {len(other_first_types)}")
print(f"Common 'first_type' between CORT and OTHER: {len(common_types)}")
if common_types:
    print("Common types:", sorted(list(common_types)))
else:
    print("No common 'first_type' cell types found.")

print(f"'first_type' unique to CORT: {len(unique_to_cort)}")
if unique_to_cort:
    print("Unique to CORT:", sorted(list(unique_to_cort)))
else:
    print("No 'first_type' unique to CORT.")

print(f"'first_type' unique to OTHER: {len(unique_to_other)}")
if unique_to_other:
    print("Unique to OTHER:", sorted(list(unique_to_other)))
else:
    print("No 'first_type' unique to OTHER.")

