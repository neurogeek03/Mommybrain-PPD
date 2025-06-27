"""
Title: Plotting EdgeR results 
Description:  Comparing the expression of a pseudobulked cell type across animals
Author:   Maria Eleni Fafouti 
Date: 23-06-2025
"""
# bash

# ========== IMPORTS ==========
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# ========== PARAMETERS ==========
project_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/EdgeR'
input_dir = os.path.join(project_dir,'edger_out')

# animals = ["BC28", "BC3", "BC9", "BC15", "BC14", "BC13"]
# results_dict = {}

# for ct in animals:
#     df = pd.read_csv(f"DE_{ct}.tsv", sep="\t", index_col=0)
#     results_dict[ct] = df
        
# ========== FUNCTIONS ==========
def plot_volcano_edger(res_df, title, ax=None, logfc_thresh=0.1, fdr_thresh=0.1):
    # Extract data
    logfc = res_df['logFC'].astype(float).values
    pval = res_df['PValue'].astype(float).values
    fdr = res_df['FDR'].astype(float).values
    genes = res_df.index.values  # Or res_df['gene'].values if gene names are in a column

    # Filter out NaNs
    valid = ~np.isnan(logfc) & ~np.isnan(pval) & ~np.isnan(fdr)
    logfc, pval, fdr, genes = logfc[valid], pval[valid], fdr[valid], genes[valid]

    # Significance criteria
    sig = (np.abs(logfc) > logfc_thresh) & (fdr < fdr_thresh)

    # Plotting
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7))

    ax.scatter(logfc, -np.log10(pval), c='gray', alpha=0.6, label='Not sig.')
    ax.scatter(logfc[sig], -np.log10(pval[sig]), c='#732B8B', alpha=0.8, label='FDR < 0.1')

    ax.axvline(-logfc_thresh, color='black', linestyle='dotted', lw=1)
    ax.axvline(logfc_thresh, color='black', linestyle='dotted', lw=1)
    ax.axhline(-np.log10(0.05), color='black', linestyle='dashed', lw=1)

    ax.set_title(title, fontsize=14)
    ax.set_xlabel("log2 fold change", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.legend()

def plot_deg_barplot(input_dir, output_path="/DE_summary_barplot.png",
                     logfc_thresh=1.0, fdr_thresh=0.05, min_genes=1,
                     sort_by="total", figsize=(12, 6), horizontal=False):
    """
    Generate a barplot of number of significantly up/downregulated genes per file.

    Parameters:
    - input_dir: str, directory with *_dge_results.tsv files
    - output_path: str, where to save the plot
    - logfc_thresh: float, log2 fold change threshold
    - fdr_thresh: float, FDR cutoff
    - min_genes: int, minimum number of DE genes to include a group
    - sort_by: str, one of "total", "upregulated", "downregulated"
    - figsize: tuple, size of the figure
    - horizontal: bool, if True, plot horizontal bars
    """
    file_paths = [os.path.join(input_dir, f) for f in os.listdir(input_dir)
                  if f.endswith("_dge_results.tsv")]
    summary = []

    for filepath in sorted(file_paths):
        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)

            logfc = df["logFC"].astype(float).values
            fdr = df["FDR"].astype(float).values

            up = np.sum((logfc > logfc_thresh) & (fdr < fdr_thresh))
            down = np.sum((logfc < -logfc_thresh) & (fdr < fdr_thresh))

            if (up + down) >= min_genes:
                title = os.path.basename(filepath).replace("_dge_results.tsv", "")
                summary.append({
                    "subclass": title,
                    "upregulated": up,
                    "downregulated": down
                })

        except Exception as e:
            print(f"Skipping {filepath} due to error: {e}")

    if not summary:
        print("No valid DE files found or no DE genes met threshold.")
        return

    df_summary = pd.DataFrame(summary)
    df_summary["total"] = df_summary["upregulated"] + df_summary["downregulated"]
    df_summary = df_summary.sort_values(sort_by)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    if horizontal:
        bar1 = ax.barh(df_summary["subclass"], df_summary["downregulated"],
               label="Downregulated (OIL > CORT)", color='blue')
        bar2 = ax.barh(df_summary["subclass"], df_summary["upregulated"],
               left=df_summary["downregulated"],
               label="Upregulated (CORT > OIL)", color='red')
        ax.set_xlabel("Number of significant genes")
    else:
        x = np.arange(len(df_summary))
        width = 0.35
        ax.bar(x - width/2, df_summary["downregulated"], width, label="Downregulated (OIL > CORT)", color='blue')
        ax.bar(x + width/2, df_summary["upregulated"], width, label="Upregulated (CORT > OIL)", color='red')
        ax.set_xticks(x)
        ax.set_xticklabels(df_summary["subclass"], rotation=90, fontsize=8)
        ax.set_ylabel("Number of significant genes")

    ax.set_title(f"DE genes per subclass (|log2FC| > {logfc_thresh}, FDR < {fdr_thresh})", fontsize=12)
    ax.legend(handles=[bar1, bar2])
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.show()


# # ========== LOOPING ==========
# # Get all .tsv files
# file_paths = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith("_dge_results.tsv")]
# n_files = len(file_paths)

# # Plotting volcanos
# ncols = 5
# nrows = int(np.ceil(n_files / ncols))
# fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 4, nrows * 4))
# axes = axes.flatten()

# for i, filepath in enumerate(sorted(file_paths)):
#     df = pd.read_csv(filepath, sep="\t", index_col=0)
#     title = os.path.basename(filepath).replace("_dge_results.tsv", "")
#     print(f"Plotting volcano for: {title}")
#     plot_volcano_edger(df, title=title, ax=axes[i])

# # Hide unused subplots
# for j in range(i + 1, len(axes)):
#     fig.delaxes(axes[j])

# plt.tight_layout()
# plt.savefig("volcano_all_celltypes.png", dpi=300)
# plt.show()

# Plotting summary barplot
plot_deg_barplot(
    input_dir=input_dir,
    output_path="volcano_plots/custom_barplot.png",
    logfc_thresh=0.1,
    fdr_thresh=0.1,
    sort_by="upregulated",  # or "downregulated" or "total"
    horizontal=True         # better for long labels
)
