import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import glob

import plotnine as p9
import liana as li
import decoupler as dc
# import omnipath as op

# Import DESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Custom functions
from functions import collapse_by_gene_symbol, print_filter_debug_info

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)
figs_dir = output_base / 'figures'
data_dir = project_path / 'data'

# =================== PARAMS ===================
cort_samples = ["BC28", "BC3", "BC9"]
oil_samples = ["BC15", "BC14", "BC13"]
sample_key = 'sample'
groupby = 'subclass_name'
condition_key = 'group'
sc.settings.figdir = figs_dir

tf_regulons = data_dir / 'RAT_TF_regulons_df.csv' #installed previously on login node
net = pd.read_csv(tf_regulons)

# # =================== INPUT ===================
# #NOTE this object contains min 10 cells per sample so were good

# slide_tags_merged = glob.glob(f'{data_dir}/*.h5ad')
# adata_path = Path(slide_tags_merged[0])
# adata = sc.read_h5ad(adata_path)
# adata.obs.loc[adata.obs[sample_key].isin(cort_samples), 'group'] = 'CORT'
# adata.obs.loc[adata.obs[sample_key].isin(oil_samples), 'group'] = 'OIL'

# # =================== CHANGE VAR INDEX ===================
# # fixing duplicates by summing
# adata = collapse_by_gene_symbol(adata, gene_symbol_col='gene_symbol')
# ########### sanity check
# duplicate_counts = adata.var.index.value_counts()
# duplicates = duplicate_counts[duplicate_counts > 1]
# print('here are the number of duplicates after switching to gene names as teh index:')
# print(duplicates.sum()) #0
# ###########

# # =========ß========== COMPUTATIONS & UMAP ===================
# # Saving count data
# adata.layers["counts"] = adata.X.copy()
# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
# sc.pl.highly_variable_genes(adata, save=True)
# sc.tl.pca(adata)
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)

# # Show pre-computed UMAP
# sc.pl.umap(adata, color=[condition_key, sample_key, 'class_name', groupby], frameon=False, ncols=2, save=True)

# adata_path = data_dir / f'object_LIANA_{adata.n_obs}_pseudobulk_gene_names_var_condition.h5ad'
# adata.write(adata_path)

# =================== IMPORT PROCESSED OBJECT ===================
adata_path = data_dir / 'object_LIANA_108123_pseudobulk_gene_names_var_condition.h5ad'
adata = sc.read_h5ad(adata_path)

# =================== PSEUDOBULK ===================
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=sample_key,
    groups_col=groupby,
    layer='counts',
    mode='sum'
)
pdata

pdata.X = np.round(pdata.X).astype(int) 

# filter samples based on number of cells and counts
dc.pp.filter_samples(pdata, min_cells = 10, min_counts=1000)
plot_path = figs_dir / 'filtered_samples.png'
dc_path_1 = str(plot_path)
dc.pl.filter_samples(pdata, groupby=[sample_key, groupby], figsize=(11, 4), save=dc_path_1)

# Save pdata.X as a CSV
ctdata_df = pd.DataFrame(ctdata.X, index=ctdata.obs_names, columns=ctdata.var_names)
ctdata_df.to_csv(output_base / 'ctdata_matrix.csv')


# =================== DE ===================
dea_results = {}
quiet = True

for cell_group in pdata.obs[groupby].unique():
    # Select cell profiles
    ctdata = pdata[pdata.obs[groupby] == cell_group].copy()

    # Add this debug statement
    if cell_group == "318 Astro-NT NN":
        print(f"DEBUG: For '{cell_group}', ctdata has {ctdata.n_obs} observations (samples) after sample filtering.")
        if ctdata.n_obs == 0:
            print(f"DEBUG: ctdata is empty for '{cell_group}'. This is likely why no genes pass the filter.")

    matrix_subset = ctdata.X[:10, :10]

    df_preview = pd.DataFrame(data=matrix_subset,index=ctdata.obs_names[:10],columns=ctdata.var_names[:10])
    print(df_preview)

    # Obtain genes that pass the edgeR-like thresholds
    # NOTE: QC thresholds might differ between cell types, consider applying them by cell type
    genes = dc.pp.filter_by_expr(ctdata,
                              group=condition_key,
                              min_count=0, # a minimum number of counts in a number of samples
                              min_total_count=0, # a minimum total number of reads across samples
                              large_n=3, 
                              min_prop=0.5
                              )
    
        # Add debug statements to inspect the 'genes' mask and resulting ctdata after filtering
    if cell_group == "318 Astro-NT NN":
        if genes is not None and genes.any():
            temp_filtered_ctdata = ctdata[:, genes].copy()
            print(f"DEBUG: For '{cell_group}', the gene filtering step returned {temp_filtered_ctdata.n_vars} genes.")
            if temp_filtered_ctdata.n_vars > 0:
                print(f"DEBUG: Preview of filtered ctdata.X for '{cell_group}' (first 5x5):")
                matrix_to_print = temp_filtered_ctdata.X
                if hasattr(matrix_to_print, 'toarray'): # Check if sparse and convert
                    matrix_to_print = matrix_to_print.toarray()
                
                df_to_print = pd.DataFrame(
                    matrix_to_print[:min(5, temp_filtered_ctdata.n_obs), :min(5, temp_filtered_ctdata.n_vars)],
                    index=temp_filtered_ctdata.obs_names[:min(5, temp_filtered_ctdata.n_obs)],
                    columns=temp_filtered_ctdata.var_names[:min(5, temp_filtered_ctdata.n_vars)]
                )
                print(df_to_print)
        else:
            print(f"DEBUG: For '{cell_group}', the gene filtering step returned 0 genes or a None mask.")


    # Filter by these genes
    if genes is None or not genes.any():
        print(f"No genes passed filter for {cell_group}, skipping.")
        continue

    ctdata = ctdata[:, genes].copy()
    
    # Build DESeq2 object
    # NOTE: this data is actually paired, so one could consider fitting the patient label as a confounder
    dds = DeseqDataSet(
        adata=ctdata,
        design_factors=condition_key,
        ref_level=[condition_key, 'OIL'], # set control as reference
        refit_cooks=True,
        quiet=quiet
    )
    
    # Compute LFCs
    dds.deseq2()
    # Contrast between stim and ctrl
    stat_res = DeseqStats(dds, contrast=[condition_key, 'CORT', 'OIL'], quiet=quiet)
    stat_res.quiet = quiet
    # Compute Wald test
    stat_res.summary()
    # Shrink LFCs
    stat_res.lfc_shrink(coeff='condition_stim_vs_ctrl') # {condition_key}_cond_vs_ref
    
    dea_results[cell_group] = stat_res.results_df

print(dea_results)