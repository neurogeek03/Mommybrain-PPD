import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def create_heatmap(df_subset, title, output_path):
    """Helper function to process and plot a specific subset of data"""
    if df_subset.empty:
        print(f"Skipping {title}: No data in this subset.")
        return

    # 1. Isolate numeric weights (ignoring first_type and is_match)
    weight_df = df_subset.select_dtypes(include=[np.number])
    
    # 2. Determine dominant type and max weight for sorting
    dominant_type = weight_df.idxmax(axis=1)
    max_weight = weight_df.max(axis=1)
    
    sort_df = pd.DataFrame({
        'dominant_type': dominant_type,
        'max_weight': max_weight
    }, index=weight_df.index).sort_values(by=['dominant_type', 'max_weight'], ascending=[True, False])
    
    # Reorder and transform
    plot_data = np.log1p(weight_df.loc[sort_df.index])

    # 3. Plotting
    height = max(10, len(plot_data.columns) * 0.35)
    plt.figure(figsize=(24, height))
    
    sns.heatmap(plot_data.T, cmap="viridis", xticklabels=False, robust=True)
    
    plt.title(f"{title} (n={len(plot_data)} cells)", fontsize=16)
    plt.xlabel("Cell IDs (Sorted by Dominant Type)")
    plt.ylabel("Cell Types")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close() # Close figure to free up memory
    print(f"Saved: {output_path}")

def main(csv_path, output_dir):
    # Load data
    print(f"Reading {csv_path}...")
    df = pd.read_csv(csv_path, index_col=0)
    
    if 'is_match' not in df.columns:
        print("Error: 'is_match' column not found in CSV. Cannot split plots.")
        sys.exit(1)

    # Split into True and False groups
    df_true = df[df['is_match'] == True]
    df_false = df[df['is_match'] == False]

    # Generate the two plots
    # We strip the extension from output_dir if it was provided as a filename
    base_name = os.path.splitext(output_dir)[0]
    
    print("Generating heatmap for Matching cells (is_match=True)...")
    create_heatmap(df_true, "RCTD Weights: Matching Cells", f"{base_name}_matching.png")
    
    print("Generating heatmap for Discordant cells (is_match=False)...")
    create_heatmap(df_false, "RCTD Weights: Discordant Cells", f"{base_name}_discordant.png")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python visualize_rctd_weights.py <input_csv> <output_base_name>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])