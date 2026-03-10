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
from functions import collapse_by_gene_symbol, print_filter_debug_info, preview_X

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

# # # # =================== INPUT ===================
# # # #NOTE this object contains min 10 cells per sample so were good

# slide_tags_merged = glob.glob(f'{data_dir}/*.h5ad')
# adata_path = Path(slide_tags_merged[0])
# print(f'Using {adata_path} as the file path ...')
# adata = sc.read_h5ad(adata_path)
# adata.obs.loc[adata.obs[sample_key].isin(cort_samples), 'group'] = 'CORT'
# adata.obs.loc[adata.obs[sample_key].isin(oil_samples), 'group'] = 'OIL'

# # # =================== CHANGE VAR INDEX ===================
# # # fixing duplicates by summing
# # adata = collapse_by_gene_symbol(adata, gene_symbol_col=adata.var.index)
# # ########### sanity check
# # duplicate_counts = adata.var.index.value_counts()
# # duplicates = duplicate_counts[duplicate_counts > 1]
# # print('here are the number of duplicates after switching to gene names as teh index:')
# # print(duplicates.sum()) #0
# # # ###########

# adata.X = adata.layers["counts"]

# #TODO check if input data needs to be norm & log transformed
# # =================== COMPUTATIONS & UMAP ===================
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

# # # =================== IMPORT PROCESSED OBJECT ===================
# # adata_path = data_dir / 'object_LIANA_108123_pseudobulk_gene_names_var_condition.h5ad'
# # adata = sc.read_h5ad(adata_path)

# # =================== PSEUDOBULK ===================
# pdata = dc.pp.pseudobulk(
#     adata,
#     sample_col=sample_key,
#     groups_col=groupby,
#     layer='counts',
#     mode='sum'
# )
# pdata

# pdata.X = np.round(pdata.X).astype(int) 

# # filter samples based on number of cells and counts
# dc.pp.filter_samples(pdata, min_cells = 10, min_counts=1000)
# plot_path = figs_dir / 'filtered_samples.png'
# dc_path_1 = str(plot_path)
# dc.pl.filter_samples(pdata, groupby=[sample_key, groupby], figsize=(11, 4), save=dc_path_1)


# # =================== DE ===================
# dea_results = {}
# quiet = True

# for cell_group in pdata.obs[groupby].unique():
#     # Select cell profiles
#     ctdata = pdata[pdata.obs[groupby] == cell_group].copy()

#     # Add this debug statement
#     if cell_group == "318 Astro-NT NN":
#         deb_path = output_base / 'ctdata_318_astro.h5ad'
#         ctdata.write(deb_path)
#         print(f"DEBUG: For '{cell_group}', ctdata has {ctdata.n_obs} observations (samples) after sample filtering.")
#         if ctdata.n_obs == 0:
#             print(f"DEBUG: ctdata is empty for '{cell_group}'. This is likely why no genes pass the filter.")

#     matrix_subset = ctdata.X[:10, :10]

#     df_preview = pd.DataFrame(data=matrix_subset,index=ctdata.obs_names[:10],columns=ctdata.var_names[:10])
#     print(df_preview)

#     # # Obtain genes that pass the edgeR-like thresholds
#     # # NOTE: QC thresholds might differ between cell types, consider applying them by cell type
#     # genes = dc.pp.filter_by_expr(ctdata,
#     #                           group=condition_key,
#     #                           min_count=1, # a minimum number of counts in a number of samples
#     #                           min_total_count=1, # a minimum total number of reads across samples
#     #                           large_n=2, 
#     #                           min_prop=0.5
#     #                           )

    
#     # Manually filter genes since library functions are failing.
#     # Keep genes that are expressed in at least 3 samples.
#     if hasattr(ctdata.X, 'toarray'):
#         expression_matrix = ctdata.X.toarray()
#     else:
#         expression_matrix = ctdata.X
#         n_samples_expressed = np.count_nonzero(expression_matrix, axis=0)
#         genes = n_samples_expressed >= 3
    
#     # ctdata_filtered_genes = sc.pp.filter_genes(ctdata, min_cells=3)

#     # Filter by these genes
#     if genes is None or not genes.any():
#         print(f"No genes passed filter for {cell_group}, skipping.")
#         continue

#     ctdata = ctdata[:, genes].copy()
    
#     # Set 'OIL' as the reference level for the 'group' condition

#     ctdata.obs[condition_key] = pd.Categorical(ctdata.obs[condition_key], categories=['OIL', 'CORT'])

#     # Build DESeq2 object using the modern formulaic design
#     # NOTE: this data is actually paired, so one could consider fitting the patient label as a confounder
#     dds = DeseqDataSet(
#         adata=ctdata,
#         design=f'~ {condition_key}',
#         refit_cooks=True,
#         quiet=quiet
#     )
    
#     # Compute LFCs
#     dds.deseq2()
#     # Contrast between stim and ctrl
#     stat_res = DeseqStats(dds, contrast=[condition_key, 'CORT', 'OIL'], quiet=quiet)
#     stat_res.quiet = quiet
#     # Compute Wald test
#     stat_res.summary()
#     # Shrink LFCs
#     stat_res.lfc_shrink(coeff='group[T.CORT]')
    
#     dea_results[cell_group] = stat_res.results_df

# print(dea_results)

# # concat results across cell types
# dea_df = pd.concat(dea_results)
# dea_df = dea_df.reset_index().rename(columns={'level_0': groupby,'level_1':'index'}).set_index('index')
# dea_df.head()

# # PyDeseq Seems to intrdoce NAs for some p-values
# # NOTE: there sometimes some NaN being introduced, best to double check that, in this case it's only for a single gene, but it might be a problem.
# len(dea_df[dea_df.isna().any(axis=1)])

# dea_path = output_base / 'first_try' / 'dea_result.csv'
# dea_df.to_csv(dea_path, index=True)

# ================ AFTER creating the DE object ================
dea_path = output_base / 'first_try' / 'dea_result.csv'
dea_df = pd.read_csv(dea_path, index_col=0)
slide_tags_merged = glob.glob(f'{data_dir}/*.h5ad')
adata_path = Path(slide_tags_merged[0])
print(f'Using {adata_path} as the file path ...')
adata = sc.read_h5ad(adata_path)

lr_res = li.multi.df_to_lr(adata,
                           dea_df=dea_df,
                           resource_name='mouseconsensus', # NOTE: uses MOUSE gene symbols!
                           expr_prop=0.1, # calculated for adata as passed - used to filter interactions
                           groupby=groupby,
                           stat_keys=['stat', 'pvalue', 'padj'],
                           use_raw=False,
                           complex_col='stat', # NOTE: we use the Wald Stat to deal with complexes
                           verbose=True,
                           return_all_lrs=False,
                           )

lr_res = lr_res.sort_values("interaction_stat", ascending=False, key=abs)
lr_res.head()
lr_res_path = output_base / f'lr_res_{stamp}.csv'
# ======================= Visualization =======================

# # Let's visualize how this looks like for all interactions  (across all cell types)
# lr_res = lr_res.sort_values("interaction_stat", ascending=False)
# hist =lr_res['interaction_stat'].hist(bins=50)
# fig = hist.get_figure()
# fig.savefig("interaction_stat_histogram.png", dpi=300)

# tileplot_fig = li.pl.tileplot(liana_res=lr_res,
#                fill = 'expr',
#                label='padj',
#                label_fun = lambda x: '*' if x < 0.05 else np.nan,
#                top_n=20,
#                orderby = 'interaction_stat',
#                orderby_ascending = False,
#                orderby_absolute = False,
#                source_title='Ligand',
#                target_title='Receptor',
#                figure_size=(16, 11)
#                )
# tileplot_fig = tileplot_fig + p9.theme(axis_text_x=p9.element_text(size=8, angle=90))
# tileplot_fig_path = output_base / 'top20_tileplot_figure.png'
# tileplot_fig.save(tileplot_fig_path, dpi=300)

# dotplot_fig = li.pl.dotplot(liana_res=lr_res,
#                      colour='interaction_stat',
#                      size='ligand_pvalue',
#                      inverse_size=True,
#                      orderby='interaction_stat',
#                      orderby_ascending=False,
#                      orderby_absolute=True,
#                      top_n=10,
#                      size_range=(0.5, 4),
#                      cmap='RdBu_r', 
#                      figure_size=(25, 10)
#                      )

# dotplot_fig = dotplot_fig + labs(size="padj (-log10 transformed)")
# dotplot_fig_path =  output_base / 'TOP10_dotplot_figure.png'
# dotplot_fig.save(dotplot_fig_path)
# # customize plot
# (
#     plot
#     + p9.theme_bw(base_size=14)
#     # fill cmap blue to red, with 0 the middle
#     + p9.scale_color_cmap('RdBu_r', midpoint=0, limits=(-10, 10))
#     # rotate x
#     + p9.theme(axis_text_x=p9.element_text(angle=90), figure_size=(11, 6))

# )

# ======================= Intracellular signaling networks =======================
# utily function to select top n interactions
def select_top_n(d, n=None):
    d = dict(sorted(d.items(), key=lambda item: abs(item[1]), reverse=True))
    return {k: v for i, (k, v) in enumerate(d.items()) if i < n}

source_label = '017 CA3 Glut'
target_label = '025 CA2-FC-IG Glut'

# NOTE: We sort by the absolute value of the interaction stat
lr_stats = lr_res[lr_res['source'].isin([source_label]) & lr_res['target'].isin([target_label])].copy()
lr_stats = lr_stats.sort_values('interaction_stat', ascending=False, key=abs)

# select receptors based on interaction stats
lr_dict = lr_stats.set_index('receptor')['interaction_stat'].to_dict()
input_scores = select_top_n(lr_dict, n=10)

# First, let's transform the DEA statistics into a DF
# we will use these to estimate deregulated TF activity
dea_wide = dea_df[[groupby, 'stat']].reset_index(names='genes').pivot(index=groupby, columns='genes', values='stat')
dea_wide = dea_wide.fillna(0)
print(dea_wide)

# Run Enrichment Analysis
estimates, pvals = dc.mt.ulm(mat=dea_wide, net=net)
estimates.T.sort_values(target_label, key=abs, ascending=False).head()

# selecting the top TFs 
tf_data = estimates.copy()
tf_dict = tf_data.loc[target_label].to_dict()
output_scores = select_top_n(tf_dict, n=5)

# obtain ppi (protein-protein interaction) network
ppis = op.interactions.OmniPath().get(genesymbols = True)

ppis['mor'] = ppis['is_stimulation'].astype(int) - ppis['is_inhibition'].astype(int)
ppis = ppis[(ppis['mor'] != 0) & (ppis['curation_effort'] >= 5) & ppis['consensus_direction']] 

input_pkn = ppis[['source_genesymbol', 'mor', 'target_genesymbol']]
input_pkn.columns = ['source', 'mor', 'target']
input_pkn.head()

# convert the PPI network into a knowledge graph
prior_graph = li.mt.build_prior_network(input_pkn, input_scores, output_scores, verbose=True)

# calculating node weights 
temp = adata[adata.obs[groupby] == target_label].copy()

node_weights = pd.DataFrame(temp.X.getnnz(axis=0) / temp.n_obs, index=temp.var_names)
node_weights = node_weights.rename(columns={0: 'props'})
node_weights = node_weights['props'].to_dict()

# ===== CORNETO =====
df_res, problem = li.mt.find_causalnet(
    prior_graph, 
    input_scores, 
    output_scores, 
    node_weights,
    # penalize (max_penalty) nodes with counts in less than 0.1 of the cells
    node_cutoff=0.1, 
    max_penalty=1,
    # the penaly of those in > 0.1 prop of cells set to:
    min_penalty=0.01,
    edge_penalty=0.1,
    verbose=False,
    max_runs=50, # NOTE that this repeats the solving either until the max runs are reached
    stable_runs=10, # or until X number of consequitive stable runs are reached (i.e. no new edges are added)
    solver='gurobi' # 'scipy' is available by default, but often results in suboptimal solutions
    )

cn.methods.carnival.visualize_network(df_res)