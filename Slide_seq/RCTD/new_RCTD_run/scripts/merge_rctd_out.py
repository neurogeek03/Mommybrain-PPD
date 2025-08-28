import pandas as pd
import glob
import os

# =============== PATHS ===============
parent_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all'

samples = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
output_dir = os.path.join(parent_dir, 'mouse_RCTD_out_merged_delta_1')
os.makedirs(output_dir, exist_ok=True)

for sample in samples:
    print(f"Processing sample: {sample}...")

    # Define output path
    output_path = os.path.join(output_dir, f"delta_1_subclass_{sample}_merged_RCTD.csv")

    # Skip if already exists
    if os.path.exists(output_path):
        print(f"⏭️  Skipping {sample}: merged file already exists at {output_path}")
    else:
        # List all RCTD result CSVs for sample 
        csv_paths = glob.glob(os.path.join(parent_dir, sample, f"{sample}_subset_*_RCTD_results.csv"))

        if not csv_paths:
            print(f"❌ No CSV files found for {sample}")
        else:
            # Merge them into one big data frame
            rctd_merged = pd.concat((pd.read_csv(csv) for csv in csv_paths), ignore_index=True)

            rctd_merged.to_csv(output_path, index=False)
            print(f"✅ Saved: {output_path}")