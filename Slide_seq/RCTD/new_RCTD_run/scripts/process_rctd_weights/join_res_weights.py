# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# import os
# import csv

# # 1. Setup Paths
# base_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/csvs_from_rctd_obj/'
# weights_path = os.path.join(base_dir, 'rctd_normalized_weights.csv')
# results_path = os.path.join(base_dir, 'rctd_summary_results.csv')
# output_path = os.path.join('/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/csvs_from_rctd_obj', 'merged_rctd_data.csv')

# # 2. Load Data
# # index_col=0 ensures the cell barcodes are used as the unique identifier (index)
# df_weights = pd.read_csv(weights_path, index_col=0)
# df_results = pd.read_csv(results_path, index_col=0)

# # Store the original weight column names to use for subsetting later
# weight_columns = df_weights.columns.tolist()

# # 3. Perform Inner Join
# # This subsets only for barcodes present in both files
# merged_df = df_weights.join(df_results, how='inner', lsuffix='_weight', rsuffix='_result')

# # 4. Filtering Step: Filter by 'spot_class'
# # Change ['singlet'] to ['singlet', 'doublet_certain'] if you want both
# target_classes = ['singlet'] 
# filtered_df = merged_df[merged_df['spot_class'].isin(target_classes)]

# print(f"Original shared cells: {len(merged_df)}")
# print(f"Cells after filtering for {target_classes}: {len(filtered_df)}")

# # 5. Remove columns from rctd_summary_results.csv
# # We keep only the columns that existed in the original weights file
# final_df = filtered_df[weight_columns]

# # 6. Save the Merged File with quoted Cell IDs
# # quoting=csv.QUOTE_NONNUMERIC ensures string IDs are in ""
# final_df.to_csv(output_path, quoting=csv.QUOTE_NONNUMERIC)

# print(f"Successfully saved cleaned weights to: {output_path}")
# print(f"Final column count: {len(final_df.columns)} (Metadata columns removed)")

import pandas as pd
import os
import csv
import numpy as np

# 1. Setup Paths
base_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/csvs_from_rctd_obj/'
weights_path = os.path.join(base_dir, 'rctd_normalized_weights.csv')
results_path = os.path.join(base_dir, 'rctd_summary_results.csv')
output_path = os.path.join(base_dir, 'merged_rctd_data.csv')

# 2. Load Data
df_weights = pd.read_csv(weights_path, index_col=0)
df_results = pd.read_csv(results_path, index_col=0)

# Capture weight column names for the calculation
weight_cols = df_weights.columns.tolist()

# 3. Join and Filter for Singlets
# We join weights with only the 'first_type' and 'spot_class' from results
merged_df = df_weights.join(df_results[['first_type', 'spot_class']], how='inner')
filtered_df = merged_df[merged_df['spot_class'] == 'singlet'].copy()

# 4. Identify the "Calculated" dominant type from weights
# idxmax finds the column name with the highest value for each row
filtered_df['calculated_max_type'] = filtered_df[weight_cols].idxmax(axis=1)

# 5. Label rows where the highest weight matches the 'first_type'
# This creates a True/False column
filtered_df['is_match'] = filtered_df['calculated_max_type'] == filtered_df['first_type']

# 6. Final Clean up
# Drop the 'spot_class' and 'calculated_max_type' as we only want weights, 
# the original first_type, and our new match label.
cols_to_keep = weight_cols + ['first_type', 'is_match']
final_df = filtered_df[cols_to_keep]

# 7. Save with quoted IDs
final_df.to_csv(output_path, quoting=csv.QUOTE_NONNUMERIC)

print(f"Saved {len(final_df)} cells.")
print(f"Consistency check: {final_df['is_match'].sum()} out of {len(final_df)} cells match.")