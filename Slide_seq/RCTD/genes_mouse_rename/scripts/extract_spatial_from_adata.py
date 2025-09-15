import os
import scanpy as sc
import pandas as pd

def extract_spatial_from_adata_files(directory, output_dir="output_spatial"):
    """
    Extracts spatial coordinates from each .h5ad file in the directory and
    saves each as a separate CSV named by the sample name (before first '_').
    """
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(directory):
        if filename.endswith(".h5ad"):
            filepath = os.path.join(directory, filename)
            print(f"Processing {filename}...")
            adata = sc.read_h5ad(filepath)

            if 'X_spatial' not in adata.obsm:
                print(f"Warning: 'X_spatial' not found in {filename}, skipping.")
                continue

            spatial_coords = adata.obsm['X_spatial']
            if spatial_coords.shape[1] != 2:
                print(f"Warning: Unexpected spatial dimension in {filename}, skipping.")
                continue

            obs_names = adata.obs_names if adata.obs_names is not None else range(spatial_coords.shape[0])
            df = pd.DataFrame(spatial_coords, columns=["x", "y"], index=obs_names)

            # Extract sample name from filename before first underscore
            sample_name = filename.split('_')[0]
            output_filepath = os.path.join(output_dir, f"{sample_name}_spatial_coords.csv")

            df.to_csv(output_filepath)
            print(f"Saved spatial coordinates to {output_filepath}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract spatial coordinates from multiple AnnData files.")
    parser.add_argument("directory", help="Directory containing .h5ad files")
    parser.add_argument("--output_dir", default="output_spatial", help="Directory to save output CSV files")

    args = parser.parse_args()
    extract_spatial_from_adata_files(args.directory, args.output_dir)
