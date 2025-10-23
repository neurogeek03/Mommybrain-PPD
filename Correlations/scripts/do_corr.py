import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

def analyze_correlations(df_dir, x_col, y_col, output_dir="correlation_plots"):
    """
    For each DataFrame in dfs, compute and plot correlation between x_col and y_col.
    Saves individual plots and a summary CSV of correlation coefficients.
    """

    os.makedirs(output_dir, exist_ok=True)
    results = []

    for fpath in df_dir.glob('*.tsv'):
        celltype = fpath.stem.split("_tags")[0]
        print(f'starting with {celltype}')

        df = pd.read_csv(fpath, sep='\t')
        data = df[[x_col, y_col]]

        # Compute Pearson correlation
        corr = data[x_col].corr(data[y_col])

        # Store result
        results.append({"dataset": celltype, "correlation": corr})

        # Plot scatter
        plt.figure(figsize=(6, 5))
        sns.regplot(data=data, x=x_col, y=y_col, scatter_kws={"alpha": 0.6})
        plt.title(f"{celltype}: log2fc per gene - Correlation = {corr:.2f}")
        plt.tight_layout()

        # Save figure
        plot_path = os.path.join(output_dir, f"{celltype}_correlation.png")
        plt.savefig(plot_path, dpi=150)
        plt.close()

    # Save correlations to CSV
    corr_df = pd.DataFrame(results)
    corr_csv = os.path.join(output_dir, "correlation_summary.csv")
    corr_df.to_csv(corr_csv, index=False)

    return corr_df

# === Run over all subdirectories ===
project_dir = Path.cwd().parents[0]
input_dir = project_dir / 'out' / 'update'

for subdir in input_dir.iterdir():
    if subdir.is_dir():
        print(f"Processing: {subdir.name}")
        summary = analyze_correlations(subdir, x_col="logFC_x", y_col="logFC_y", output_dir=f'indiv_plots_{subdir.name}')

        if summary.empty:
            continue

        plt.figure(figsize=(10, 5))
        sns.barplot(
            data=summary.sort_values("correlation", ascending=False),
            x="dataset",
            y="correlation"
        )
        plt.xticks(rotation=90)
        plt.title(f"{subdir.name}: correlation coefficients")
        plt.tight_layout()
        fig_path = input_dir / f'{subdir.name}_vs_tags_summary.png'
        plt.savefig(fig_path, dpi=150)
        plt.close()
        print(f"✅ Saved plot: {fig_path}")