import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata
import logging
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Constants ---
CELL_TYPE_COL = 'subclass'
REGION_COL = 'region'

def stratified_subsample(adata, n_obs, groupby, min_per_group, random_state=0):
    """
    Subsamples an AnnData object to a target number of observations, ensuring a
    minimum number of observations per group are retained.
    """
    np.random.seed(random_state)

    if n_obs >= adata.n_obs:
        return adata.copy()

    obs_df = adata.obs.copy()
    
    # Determine which cells to keep to satisfy the minimum constraint
    protected_indices = []
    for group_name, group_df in obs_df.groupby(groupby):
        if len(group_df) < min_per_group:
            protected_indices.extend(group_df.index)
        else:
            protected_indices.extend(
                np.random.choice(group_df.index, size=min_per_group, replace=False)
            )
    
    protected_indices = list(set(protected_indices))
    n_protected = len(protected_indices)
    
    if n_protected >= n_obs:
        logging.warning(f"The number of cells to meet minimums ({n_protected}) "
                        f"exceeds the target ({n_obs}). Downsampling will not guarantee all minimums. "
                        "Sampling {n_obs} cells from the protected set.")
        final_indices = np.random.choice(protected_indices, size=n_obs, replace=False)
    else:
        final_indices = protected_indices
        remaining_indices = list(set(obs_df.index) - set(protected_indices))
        n_to_sample_more = n_obs - n_protected
        
        if n_to_sample_more > 0 and len(remaining_indices) > 0:
            more_indices = np.random.choice(remaining_indices, size=n_to_sample_more, replace=False)
            final_indices.extend(more_indices)

    return adata[final_indices, :].copy()


def process_region(filepath, meta_df, qc_plot_dir, min_genes, min_counts, target_cells, min_cells_per_type):
    """Loads, plots QC, filters, merges metadata, and downsamples a single region."""
    try:
        filename = os.path.basename(filepath)
        brain_region = filename.split('-', 2)[-1].replace("-raw.h5ad", "")
        logging.info(f"🧠 Processing region: {brain_region}")

        ad = sc.read_h5ad(filepath)
        logging.info(f"Loaded {ad.n_obs} cells and {ad.n_vars} genes.")

        # --- Plot pre-filter QC metrics ---
        sc.pp.calculate_qc_metrics(ad, percent_top=None, log1p=False, inplace=True)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        sns.violinplot(y=ad.obs['total_counts'], ax=axes[0]).set_title('Counts per Cell')
        sns.violinplot(y=ad.obs['n_genes_by_counts'], ax=axes[1]).set_title('Genes per Cell')
        fig.suptitle(f'QC Metrics for {brain_region} (before filtering)')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plot_path = os.path.join(qc_plot_dir, f"{brain_region}_qc_violin.png")
        plt.savefig(plot_path, dpi=150)
        plt.close(fig)
        logging.info(f"Saved QC plot to {plot_path}")

        # --- Filter low-quality cells ---
        sc.pp.filter_cells(ad, min_genes=min_genes)
        sc.pp.filter_cells(ad, min_counts=min_counts)
        if ad.n_obs == 0:
            logging.warning(f"No cells remaining for {brain_region} after filtering. Skipping.")
            return None
        logging.info(f"Filtered to {ad.n_obs} cells.")

        # --- Merge Metadata (before downsampling) ---
        ad.obs[REGION_COL] = brain_region
        ad.obs.index.name = 'cell_label'
        meta_region = meta_df[meta_df[REGION_COL] == brain_region].set_index('cell_label')
        new_meta_cols = meta_region.columns.difference(ad.obs.columns)
        ad.obs = ad.obs.join(meta_region[new_meta_cols], how='left')
        
        if CELL_TYPE_COL not in ad.obs.columns:
            logging.warning(f"Cell type column '{CELL_TYPE_COL}' not found for region {brain_region}. Skipping.")
            return None
        ad = ad[ad.obs[CELL_TYPE_COL].notna()].copy()
        logging.info(f"{ad.n_obs} cells remain after removing cells with no subclass annotation.")

        # --- Stratified Downsampling ---
        if ad.n_obs > target_cells:
            logging.info(f"Region has {ad.n_obs} cells, performing stratified downsampling to {target_cells}.")
            ad = stratified_subsample(ad, target_cells, CELL_TYPE_COL, min_cells_per_type, random_state=0)
            logging.info(f"Downsampled to {ad.n_obs} cells.")

        return ad

    except Exception as e:
        logging.error(f"Failed to process {filepath}: {e}", exc_info=True)
        return None


def main():
    parser = argparse.ArgumentParser(description="Merge, filter, and annotate regional cell type data.")
    parser.add_argument('--in_dir', type=str, default='/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/data', help='Input directory with .h5ad files.')
    parser.add_argument('--out_dir', type=str, default='/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference_coronal', help='Output directory.')
    parser.add_argument('--metadata_csv', type=str, default='/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/subset_cellids_celltypes_metadata_20241115.csv', help='Path to metadata CSV.')
    parser.add_argument('--min_genes', type=int, default=200, help='Minimum genes per cell.')
    parser.add_argument('--min_counts', type=int, default=500, help='Minimum counts per cell.')
    parser.add_argument('--target_cells_per_region', type=int, default=50000, help='Max cells per region after sampling.')
    parser.add_argument('--min_cells_per_type', type=int, default=25, help='Minimum cells per type to preserve during sampling.')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    qc_plot_dir = os.path.join(args.out_dir, 'qc_plots')
    os.makedirs(qc_plot_dir, exist_ok=True)

    # --- Load and Prepare Metadata ---
    logging.info("Loading and preparing metadata...")
    meta_df = pd.read_csv(args.metadata_csv)
    
    # --- Conditional Region Parsing ---
    # The logic for extracting the region name from 'feature_matrix_label'
    # differs for Isocortex vs. other regions.
    logging.info("Parsing region names from 'feature_matrix_label'...")

    # Create a boolean mask for rows related to Isocortex
    isocortex_mask = meta_df["feature_matrix_label"].str.contains("Isocortex", na=False)

    # For Isocortex rows: split into 3 parts, take the last
    isocortex_regions = meta_df.loc[isocortex_mask, "feature_matrix_label"].str.split("-", n=2).str[2]

    # For other rows: split on all '-', take the 3rd part
    other_regions = meta_df.loc[~isocortex_mask, "feature_matrix_label"].str.split("-").str[2]

    # Combine the results back into the 'region' column and sort by original index
    meta_df[REGION_COL] = pd.concat([isocortex_regions, other_regions]).sort_index()

    logging.info("Example of parsed regions from metadata:")
    logging.info(f"\n{meta_df[REGION_COL].value_counts().head().to_string()}")
    logging.info(f"Metadata loaded for {meta_df.shape[0]} cells.")

    logging.info(f"Metadata loaded for {meta_df.shape[0]} cells.")

    # --- Process Each Region ---
    adata_list = []
    for filename in sorted(os.listdir(args.in_dir)):
        if filename.endswith('.h5ad'):
            filepath = os.path.join(args.in_dir, filename)
            processed_ad = process_region(
                filepath, meta_df, qc_plot_dir, args.min_genes, args.min_counts,
                args.target_cells_per_region, args.min_cells_per_type
            )
            if processed_ad:
                adata_list.append(processed_ad)

    if not adata_list:
        logging.error("No AnnData objects were processed successfully. Exiting.")
        return

    # --- Concatenate and Final Steps ---
    logging.info("🔗 Concatenating all regions...")
    adata_all = anndata.concat(adata_list, join='outer', merge='same', index_unique=None)
    logging.info(f"Final concatenated shape: {adata_all.shape}")

    # --- Final filtering on the combined object (as a sanity check) ---
    logging.info(f"🔬 Final filtering of cell types in '{CELL_TYPE_COL}' with less than {args.min_cells_per_type} cells globally...")
    cell_counts = adata_all.obs[CELL_TYPE_COL].value_counts()
    valid_cell_types = cell_counts[cell_counts >= args.min_cells_per_type].index
    n_obs_before = adata_all.n_obs
    adata_all = adata_all[adata_all.obs[CELL_TYPE_COL].isin(valid_cell_types)].copy()
    logging.info(f"Filtered from {n_obs_before} to {adata_all.n_obs} cells globally.")

    # --- Create and Save Summary CSV ---
    logging.info("📄 Creating and saving summary CSV...")
    summary_df = adata_all.obs.reset_index()[['cell_label', REGION_COL, CELL_TYPE_COL]]
    subclass_counts = summary_df[CELL_TYPE_COL].map(summary_df[CELL_TYPE_COL].value_counts())
    summary_df['total_cells_per_subclass'] = subclass_counts
    summary_df.rename(columns={CELL_TYPE_COL: 'subclass'}, inplace=True)
    csv_output_path = os.path.join(args.out_dir, "cell_subclass_summary.csv")
    summary_df.to_csv(csv_output_path, index=False)
    logging.info(f"✅ Saved summary to: {csv_output_path}")

    # --- Save Final AnnData Object ---
    output_filename = "celltypes_WMB_high_qual_cells_all_regions_max50k.h5ad"
    output_filepath = os.path.join(args.out_dir, output_filename)
    adata_all.write_h5ad(output_filepath)
    logging.info(f"✅ Saved final AnnData to: {output_filepath}")

if __name__ == '__main__':
    main()