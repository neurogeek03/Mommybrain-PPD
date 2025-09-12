import pandas as pd
import argparse 
import glob
import os

# =============== ARGS ===============
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_dir", required=True, help="Path to output folder to store the merged matadata csvs")
args = parser.parse_args()

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

print(f"Results will be written to: {output_dir}")

# =============== SCRIPT ===============
parent_dir = '/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all'

samples = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
# samples =["B37"]

os.makedirs(output_dir, exist_ok=True)

for sample in samples:
    print(f"Processing sample: {sample}...")

    # Define output path
    output_path = os.path.join(output_dir, f"delta_5_umi30_{sample}_merged_RCTD.csv")

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
            print(f"Saved: {output_path}")