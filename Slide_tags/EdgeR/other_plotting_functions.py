# MORE SPECIFIC FUNCTIONS
def plot_and_save_individual_volcanos(input_dir, output_dir, logfc_thresh=0.1, fdr_thresh=0.1, top_n=0):
    """
    Creates and saves individual volcano plots for each DE result file, labeling subclass-specific genes only.

    Parameters:
    - input_dir: str, directory containing *_dge_results.tsv files
    - output_dir: str, directory where individual plots will be saved
    - logfc_thresh: float, log2 fold change threshold for significance
    - fdr_thresh: float, FDR threshold for significance
    - top_n: int, ignored (labels are defined by subclass)
    """
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    # Gene labels for specific subclasses
    subclass_genes = {
        "327_Oligo_NN": ["Fkbp5", "Slc2a1", "Nfkbia", "Bcl2l1", "Aspa"],
        "007_L2_3_IT_CTX_Glut": ["Fkbp5", "Btg2", "Camk1g", "Zfp189", "Il6r"],
        "061_STR_D1_Gaba": ["Camk1g", "Zfp189", "Actg1", "Baalc", "Tmem196"]
    }

    os.makedirs(output_dir, exist_ok=True)

    file_paths = [os.path.join(input_dir, f)
                  for f in os.listdir(input_dir)
                  if f.endswith("_dge_results.tsv")]

    for filepath in sorted(file_paths):
        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)
            df = df.dropna(subset=["logFC", "PValue", "FDR"])

            title = os.path.basename(filepath).replace("_dge_results.tsv", "")
            fig, ax = plt.subplots(figsize=(8, 6))

            # Plot volcano
            plot_volcano_edger(df, title=title, ax=ax,
                               logfc_thresh=logfc_thresh,
                               fdr_thresh=fdr_thresh)

            # Label only predefined genes for this subclass
            if title in subclass_genes:
                for gene in subclass_genes[title]:
                    if gene in df.index:
                        row = df.loc[gene]
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


def plot_and_save_individual_volcanos(input_dir, output_dir, logfc_thresh=0.1, fdr_thresh=0.1, top_n=0):
    """
    Creates and saves individual volcano plots for each DE result file, labeling subclass-specific genes only.

    Parameters:
    - input_dir: str, directory containing *_dge_results.tsv files
    - output_dir: str, directory where individual plots will be saved
    - logfc_thresh: float, log2 fold change threshold for significance
    - fdr_thresh: float, FDR threshold for significance
    - top_n: int, ignored (labels are defined by subclass)
    """

    # Gene labels for specific subclasses
    subclass_genes = {
        "327_Oligo_NN": ["Fkbp5", "Slc2a1", "Nfkbia", "Bcl2l1", "Aspa"],
        "007_L2_3_IT_CTX_Glut": ["Fkbp5", "Btg2", "Camk1g", "Zfp189", "Il6r"],
        "061_STR_D1_Gaba": ["Camk1g", "Zfp189", "Actg1", "Baalc", "Tmem196"]
    }

    os.makedirs(output_dir, exist_ok=True)

    file_paths = [os.path.join(input_dir, f)
                  for f in os.listdir(input_dir)
                  if f.endswith("_dge_results.tsv")]

    for filepath in sorted(file_paths):
        try:
            df = pd.read_csv(filepath, sep="\t", index_col=0)
            df = df.dropna(subset=["logFC", "PValue", "FDR"])

            title = os.path.basename(filepath).replace("_dge_results.tsv", "")
            fig, ax = plt.subplots(figsize=(8, 6))

            # Plot volcano
            plot_volcano_edger(df, title=title, ax=ax,
                               logfc_thresh=logfc_thresh,
                               fdr_thresh=fdr_thresh)

            # Label only predefined genes for this subclass
            if title in subclass_genes:
                for gene in subclass_genes[title]:
                    if gene in df.index:
                        row = df.loc[gene]
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