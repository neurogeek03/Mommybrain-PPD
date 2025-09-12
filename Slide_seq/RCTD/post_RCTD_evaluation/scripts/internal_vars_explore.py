import pandas as pd
from tqdm import tqdm   # <--- Make sure this is imported

df = pd.read_csv("/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/post_RCTD_evaluation/RCTD_object_QC/internal_vars_summary_wide.csv")
df.head(10)

# Filter to one file / subset
subset = df[df['file'] == 'B03_subset_1_RCTD.rds'].iloc[0]  # first row of this subset

print(subset)

cols_to_check = [col for col in df.columns if col not in ["sigma", "file", "cell_types_assigned"]]

# # Convert string columns to lists and compute lengths
# lengths_df = pd.DataFrame()
# for col in cols_to_check:
#     # Convert to list (even if only one element)
#     lengths_df[f"{col}_length"] = df[col].astype(str).str.split(",").apply(len)

# # Keep metadata
# lengths_df["cell_types_assigned"] = df["cell_types_assigned"]
# lengths_df["sigma"] = df["sigma"]
# lengths_df["file"] = df["file"]

lengths = {}
for col in cols_to_check:
    # Convert to list and get length
    value = subset[col]
    if pd.isna(value):
        lengths[f"{col}_length"] = 0
    else:
        lengths[f"{col}_length"] = len(str(value).split(","))

# Add metadata
lengths["cell_types_assigned"] = subset["cell_types_assigned"]
lengths["sigma"] = subset["sigma"]
lengths["file"] = subset["file"]

# Convert to DataFrame for easier viewing / saving
lengths_df = pd.DataFrame([lengths])

# Preview
pd.set_option('display.max_columns', None)
print(lengths_df.head())

# Save CSV
df.to_csv("B03_subset_1_internal_vars_unpacked_progress.csv", index=False)