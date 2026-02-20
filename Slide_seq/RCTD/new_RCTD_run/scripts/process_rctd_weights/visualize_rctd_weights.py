import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def main(csv_path, output_path):
    # Load data
    print(f"Reading {csv_path}...")
    df = pd.read_csv(csv_path, index_col=0)
    
    # Determine the dominant cell type for each barcode to use for sorting
    # This helps group similar cells together on the X-axis
    print("Ordering cells by dominant cell type...")
    dominant_type = df.idxmax(axis=1)
    max_weight = df.max(axis=1)
    
    # Create a temporary sorting dataframe
    sort_df = pd.DataFrame({
        'dominant_type': dominant_type,
        'max_weight': max_weight
    }, index=df.index)
    
    # Sort by type name, then by weight magnitude within that type
    sort_df = sort_df.sort_values(by=['dominant_type', 'max_weight'], ascending=[True, False])
    
    # Reorder the original dataframe
    df_sorted = df.loc[sort_df.index]
    
    # Use log1p values for the weights to handle zeros and small values
    print("Applying log1p transformation...")
    df_log = np.log1p(df_sorted)
    
    print(f"Data shape: {df_log.shape}")

    # To make the X-axis (index) legible, we usually don't show all labels.
    # Adjust figure height based on number of cell types to keep Y-axis legible
    height = max(10, len(df_log.columns) * 0.15)
    plt.figure(figsize=(24, height))
    
    # Plotting Transposed so Cell Types are on Y-axis and Cells are on X-axis
    ax = sns.heatmap(df_log.T, cmap="viridis", xticklabels=False)
    
    plt.title("RCTD Cell Type Weights (log1p scale, grouped by dominant type)")
    plt.xlabel("Cell IDs (Barcodes sorted by dominant cell type)")
    plt.ylabel("Cell Types")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved heatmap to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python visualize_rctd_weights.py <input_csv> <output_png>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])
