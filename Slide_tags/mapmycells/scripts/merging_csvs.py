import pandas as pd
import os
import glob

# Specify the directory containing the CSV files
csv_dir = "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/Post_bender"

sample_list = ['BC9', 'BC15', 'BC14', 'BC28']

for sample in sample_list: 
    print(f"Processing {sample}...")

    sample_folder = os.path.join(csv_dir, sample)

    # Find all CSV files in the directory
    csv_files = glob.glob(os.path.join(sample_folder, "*.csv"))

    print(f"Found {len(csv_files)} CSV files.")

    # Read and merge CSVs
    dataframes = []
    for csv_path in csv_files:
        try:
            df = pd.read_csv(csv_path, skiprows=4, index_col=0)
            dataframes.append(df)
        except Exception as e:
            print(f"Error reading {csv_path}: {e}")

    # Concatenate all dataframes
    merged_df = pd.concat(dataframes, axis=0)

    print(f"Merged DataFrame shape: {merged_df.shape}")

    output_path = os.path.join(sample_folder, f"{sample}_mapmycells.csv")

    merged_df.to_csv(output_path)

    print(f"Merged CSV saved to: {output_path}")
