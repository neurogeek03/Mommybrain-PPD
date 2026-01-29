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
from functions import collapse_by_gene_symbol

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

# =================== INPUT ===================
#TODO check filtering criteria for DE sample 

#TODO object needs to be collapsed to hgene name 
slide_tags_merged = glob.glob(f'{data_dir}/*.h5ad')
adata_path = Path(slide_tags_merged[0])
adata = sc.read_h5ad(adata_path)

# =================== CHANGE VAR INDEX ===================
# adata.var = adata.var.set_index('gene_symbol')
# adata.var.index = adata.var.index.astype(str)

adata = collapse_by_gene_symbol(adata, gene_symbol_col='gene_symbol')
# sanity check
duplicate_counts = adata.var.index.value_counts()
duplicates = duplicate_counts[duplicate_counts > 1]
print('here are the number of duplicates after switching to gene names as teh index:')
print(duplicates.sum()) #242
####

#TODO might need to use collapse function
# fixing duplicates by summing
# adata = aggregate_duplicate_genes(adata)


# sanity check
duplicate_counts = adata.var.index.value_counts()
duplicates = duplicate_counts[duplicate_counts > 1]
print('here are the number of duplicates after switching to gene names as teh index:')
print(duplicates.sum()) #242
####


tf_regulons = data_dir / 'RAT_TF_regulons_df.csv' #installed previously on login node
net = pd.read_csv(tf_regulons)

adata.obs.loc[adata.obs[sample_key].isin(cort_samples), 'group'] = 'CORT'
adata.obs.loc[adata.obs[sample_key].isin(oil_samples), 'group'] = 'OIL'

# computations
# # Saving count data
# adata.layers["counts"] = adata.X.copy()
# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
# sc.pl.highly_variable_genes(adata, save=True)
# sc.tl.pca(adata)
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)


# # =================== SHOWCASE ===================
# # Show pre-computed UMAP
# sc.pl.umap(adata, color=[condition_key, sample_key, 'class_name', groupby], frameon=False, ncols=2, save=True)

# =================== PSEUDOBULK ===================
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=sample_key,
    groups_col=groupby,
    layer='counts',
    mode='sum'
)
pdata

# filter samples based on number of cells and counts
dc.pp.filter_samples(pdata, min_cells = 10, min_counts=1000)
plot_path = figs_dir / 'filtered_samples.png'
dc_path_1 = str(plot_path)
dc.pl.filter_samples(pdata, groupby=[sample_key, groupby], figsize=(11, 4), save=dc_path_1)


# =================== DE ===================
dea_results = {}
quiet = True

for cell_group in pdata.obs[groupby].unique():
    # Select cell profiles
    ctdata = pdata[pdata.obs[groupby] == cell_group].copy()

    # Obtain genes that pass the edgeR-like thresholds
    # NOTE: QC thresholds might differ between cell types, consider applying them by cell type
    genes = dc.pp.filter_by_expr(ctdata,
                              group=condition_key,
                              min_count=5, # a minimum number of counts in a number of samples
                              min_total_count=10 # a minimum total number of reads across samples
                              )

    # Filter by these genes
    ctdata = ctdata[:, genes].copy()
    
    # Build DESeq2 object
    # NOTE: this data is actually paired, so one could consider fitting the patient label as a confounder
    dds = DeseqDataSet(
        adata=ctdata,
        design_factors=condition_key,
        ref_level=[condition_key, 'ctrl'], # set control as reference
        refit_cooks=True,
        quiet=quiet
    )
    
    # Compute LFCs
    dds.deseq2()
    # Contrast between stim and ctrl
    stat_res = DeseqStats(dds, contrast=[condition_key, 'stim', 'ctrl'], quiet=quiet)
    stat_res.quiet = quiet
    # Compute Wald test
    stat_res.summary()
    # Shrink LFCs
    stat_res.lfc_shrink(coeff='condition_stim_vs_ctrl') # {condition_key}_cond_vs_ref
    
    dea_results[cell_group] = stat_res.results_df

print(dea)