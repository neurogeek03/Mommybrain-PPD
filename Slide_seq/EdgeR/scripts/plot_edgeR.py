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
# project_dir = '/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/EdgeR'
# input_dir = os.path.join(project_dir,'example_plots')
# output_dir = os.path.join(project_dir, 'figures')

# animals = ["BC28", "BC3", "BC9", "BC15", "BC14", "BC13"]
# results_dict = {}

# for ct in animals:
#     df = pd.read_csv(f"DE_{ct}.tsv", sep="\t", index_col=0)
#     results_dict[ct] = df
        
# ========== FUNCTIONS ==========
def plot_volcano_edger(res_df, title, ax=None, logfc_thresh=0.1, fdr_thresh=0.1):
    """
    Plots a volcano plot for differential expression analysis results.

    Parameters:
    res_df (DataFrame): A pandas DataFrame containing the results of the differential expression analysis.
                        It must contain 'logFC', 'PValue', and 'FDR' columns.
    title (str): The title of the plot.
    ax (matplotlib.axes.Axes, optional): The axes on which to plot. If None, a new figure and axes will be created.
    logfc_thresh (float, optional): The log2 fold change threshold for significance. Default is 0.1.
    fdr_thresh (float, optional): The false discovery rate threshold for significance. Default is 0.1.

    Returns:
    None: The function displays the plot.
    """
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
                     sort_by="total", figsize=(12, 8), horizontal=False,
                     log_scale=False):
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
    - log_scale: bool, if True, use log scale for the count axis
    """
    file_paths = [os.path.join(input_dir, f) for f in os.listdir(input_dir)
                  if f.endswith("_edgeR_results.tsv")]
    summary = []

    for filepath in sorted(file_paths):
        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)

            logfc = df["logFC"].astype(float).values
            fdr = df["FDR"].astype(float).values

            up = np.sum((logfc > logfc_thresh) & (fdr < fdr_thresh))
            down = np.sum((logfc < -logfc_thresh) & (fdr < fdr_thresh))

            if (up + down) >= min_genes:
                title = os.path.basename(filepath).replace("_edgeR_results.tsv", "")
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

    # --- Get number of bars for margin calculation ---
    num_bars = len(df_summary)
    if num_bars == 0:
        print("No subclasses to plot after filtering.")
        return

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    if horizontal:
        bar1 = ax.barh(df_summary["subclass"], df_summary["downregulated"],
               label="Downregulated (OIL > CORT)", color='blue')
        bar2 = ax.barh(df_summary["subclass"], df_summary["upregulated"],
               left=df_summary["downregulated"],
               label="Upregulated (CORT > OIL)", color='red')
        ax.set_xlabel("Number of significant genes")
        if log_scale:
            ax.set_xscale("log")
        
        # --- START: Added code for HORIZONTAL plot ---
        # 1. Remove top and bottom whitespace
        ax.set_ylim(num_bars - 0.5, -0.5)
        
        # 2. Remove plot borders (spines), keep x-axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False) # Hide the y-axis line
        ax.tick_params(axis='y', length=0) # Hide the y-axis ticks
        # --- END: Added code ---

    else:
        x = np.arange(num_bars)
        width = 0.35
        bar1 = ax.bar(x - width/2, df_summary["downregulated"], width, label="Downregulated (OIL > CORT)", color='blue')
        bar2 = ax.bar(x + width/2, df_summary["upregulated"], width, label="Upregulated (CORT > OIL)", color='red')
        ax.set_xticks(x)
        ax.set_xticklabels(df_summary["subclass"], rotation=90, fontsize=8)
        ax.set_ylabel("Number of significant genes")
        if log_scale:
            ax.set_yscale("log")
            
        # --- START: Added code for VERTICAL plot ---
        # 1. Remove left and right whitespace
        ax.set_xlim(-0.5, num_bars - 0.5)

        # 2. Remove plot borders (spines), keep y-axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False) # Hide the x-axis line
        ax.tick_params(axis='x', length=0) # Hide the x-axis ticks
        # --- END: Added code ---
    # # Plot
    # fig, ax = plt.subplots(figsize=figsize)
    # if horizontal:
    #     bar1 = ax.barh(df_summary["subclass"], df_summary["downregulated"],
    #            label="Downregulated (OIL > CORT)", color='blue')
    #     bar2 = ax.barh(df_summary["subclass"], df_summary["upregulated"],
    #            left=df_summary["downregulated"],
    #            label="Upregulated (CORT > OIL)", color='red')
    #     ax.set_xlabel("Number of significant genes")
    #     if log_scale:
    #         ax.set_xscale("log")
    # else:
    #     x = np.arange(len(df_summary))
    #     width = 0.35
    #     bar1 = ax.bar(x - width/2, df_summary["downregulated"], width, label="Downregulated (OIL > CORT)", color='blue')
    #     bar2 = ax.bar(x + width/2, df_summary["upregulated"], width, label="Upregulated (CORT > OIL)", color='red')
    #     ax.set_xticks(x)
    #     ax.set_xticklabels(df_summary["subclass"], rotation=90, fontsize=8)
    #     ax.set_ylabel("Number of significant genes")
    
    ax.set_title(f"DE genes per subclass (|log2FC| > {logfc_thresh}, FDR < {fdr_thresh})", fontsize=12)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)

def plot_and_save_individual_volcanos(input_dir, output_dir, logfc_thresh=0.1, fdr_thresh=0.1, top_n=10):
    """
    Creates and saves individual volcano plots for each DE result file, labeling top N genes by FDR.

    Parameters:
    - input_dir: str, directory containing *_dge_results.tsv files
    - output_dir: str, directory where individual plots will be saved
    - logfc_thresh: float, log2 fold change threshold for significance
    - fdr_thresh: float, FDR threshold for significance
    - top_n: int, number of top genes (lowest FDR) to label
    """
    os.makedirs(output_dir, exist_ok=True)

    file_paths = [os.path.join(input_dir, f)
                  for f in os.listdir(input_dir)
                  if f.endswith("_edgeR_results.tsv")]

    for filepath in sorted(file_paths):
        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)
            df = df.dropna(subset=["logFC", "PValue", "FDR"])

            title = os.path.basename(filepath).replace("_dge_results.tsv", "")

            # Sort by FDR for top genes
            top_genes = df.nsmallest(top_n, "FDR")

            fig, ax = plt.subplots(figsize=(8, 6))
            plot_volcano_edger(df, title=title, ax=ax, logfc_thresh=logfc_thresh, fdr_thresh=fdr_thresh)

            # Label top genes
            for gene, row in top_genes.iterrows():
                ax.annotate(
                    gene,
                    (row["logFC"], -np.log10(row["PValue"])),
                    fontsize=7,
                    xytext=(4, 4),
                    textcoords="offset points",
                    arrowprops=dict(arrowstyle='-', lw=0.5),
                )

            plt.tight_layout()
            outpath = os.path.join(output_dir, f"{title}_volcano.png")
            plt.savefig(outpath, dpi=300)
            plt.close()
            print(f"Saved: {outpath}")

        except Exception as e:
            print(f"Skipping {filepath} due to error: {e}")




# # ========== COMBINED VOLCANO PLOTS ==========
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
# plt.savefig("filtered_volcano_all_celltypes.png", dpi=300)
# plt.show()

# # ========== INDIVIDUAL VOLCANO PLOTS ==========
# plot_and_save_individual_volcanos(
#     input_dir=input_dir,
#     output_dir="figures/nolabels_individual_filtered",
#     logfc_thresh=0.1,
#     fdr_thresh=0.1,
#     top_n=5
# )


# # ========== SUMMARY BARPLOT ==========
# plot_deg_barplot(
#     input_dir=input_dir,
#     output_path="volcano_plots/filtered_custom_barplot.png",
#     logfc_thresh=0.1,
#     fdr_thresh=0.1,
#     sort_by="upregulated",  # or "downregulated" or "total"
#     horizontal=True         # better for long labels
# )

# def plot_and_save_individual_volcanos(input_dir, output_dir, logfc_thresh=0.1, fdr_thresh=0.1, top_n=0):
#     """
#     Creates and saves individual volcano plots for each DE result file, labeling subclass-specific genes only.

#     Parameters:
#     - input_dir: str, directory containing *_dge_results.tsv files
#     - output_dir: str, directory where individual plots will be saved
#     - logfc_thresh: float, log2 fold change threshold for significance
#     - fdr_thresh: float, FDR threshold for significance
#     - top_n: int, ignored (labels are defined by subclass)
#     """

#     # Gene labels for specific subclasses
#     subclass_genes = {
#     "327_Oligo_NN": ["Fkbp5", "Sgk1", "Ddit4", "Pdk4", "Btg2"],
#     "007_L2_3_IT_CTX_Glut": ["Fkbp5", "Btg2", "Zfp189", "Tgfb2", "Wipf3"],
#     "061_STR_D1_Gaba": ["Camk1g", "Actg1", "Kcnh7", "Arhgap26"]
#     }


#     os.makedirs(output_dir, exist_ok=True)

#     file_paths = [os.path.join(input_dir, f)
#                   for f in os.listdir(input_dir)
#                   if f.endswith("_dge_results.tsv")]

#     for filepath in sorted(file_paths):
#         try:
#             df = pd.read_csv(filepath, sep="\t", index_col=0)
#             df = df.dropna(subset=["logFC", "PValue", "FDR"])

#             title = os.path.basename(filepath).replace("_dge_results.tsv", "")
#             fig, ax = plt.subplots(figsize=(8, 6))

#             # Plot volcano
#             plot_volcano_edger(df, title=title, ax=ax,
#                                logfc_thresh=logfc_thresh,
#                                fdr_thresh=fdr_thresh)

#             # Label only predefined genes for this subclass
#             if title in subclass_genes:
#                 for gene in subclass_genes[title]:
#                     if gene in df.index:
#                         row = df.loc[gene]
#                         ax.annotate(
#                             gene,
#                             (row["logFC"], -np.log10(row["PValue"])),
#                             fontsize=10,
#                             xytext=(4, 4),
#                             textcoords="offset points",
#                             arrowprops=dict(arrowstyle='-', lw=0.5),
#                         )

#             plt.tight_layout()
#             outpath = os.path.join(output_dir, f"{title}_volcano.png")
#             plt.savefig(outpath, dpi=300)
#             plt.close()
#             print(f"Saved: {outpath}")

#         except Exception as e:
#             print(f"Skipping {filepath} due to error: {e}")


if __name__ == "__main__":
    # Path to the main results directory
    base_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR/out/edger_out"
    # List all subfolders (comparisons)
    comparison_folders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]
    print(comparison_folders)
    for comp in comparison_folders:
        comp_dir = os.path.join(base_dir, comp)
        # output_dir = os.path.join(base_dir, 'new_01_fdr_threshold', comp)

        output_path = os.path.join(base_dir, '1M_cells', f"{comp}_barplot.png")
        print(f"Plotting barplot for comparison: {comp}")
        plot_deg_barplot(
            input_dir=comp_dir,
            output_path=output_path,
            logfc_thresh=0.1,   # accept any fold change
            fdr_thresh=0.1,     # accept any FDR
            sort_by="total",
            figsize=(10, 8),
            # log_scale = True,
            horizontal=True
        )
        print(f"saved results to {output_path} ")

    # base_dir = "/scratch/mfafouti/Mommybrain/Slide_seq/EdgeR/edger_out/02_fdr_threshold"
    # comparison_folders = [
    #     f for f in os.listdir(base_dir)
    #     if os.path.isdir(os.path.join(base_dir, f))
    # ]

    # for comp in comparison_folders:
    #     comp_dir = os.path.join(base_dir, comp)

    #     # Each comparison gets its own output folder under "plots_volcano"
    #     output_dir = os.path.join(base_dir, "plots_volcano", comp)
    #     os.makedirs(output_dir, exist_ok=True)

    #     print(f"Plotting volcano plots for comparison: {comp}")
    #     plot_and_save_individual_volcanos(
    #         input_dir=comp_dir,
    #         output_dir=output_dir,   # <-- directory, not file
    #         logfc_thresh=0.1,
    #         fdr_thresh=0.1,
    #         top_n=30
    #     )

