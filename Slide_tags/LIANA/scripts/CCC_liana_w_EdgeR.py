from datetime import datetime
from pathlib import Path
import pandas as pd
import scanpy as sc
import glob

import liana as li
# import omnipath as op

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
print(f"current working directory: {project_path}")
output_base = project_path / 'out'
output_base.mkdir(exist_ok=True, parents=True)
figs_dir = output_base / 'figures'
data_dir = project_path / 'data'

dea_path = output_base / 'edgeR_dge' / 'edgeR_dge_input_liana.csv'

# =================== PARAMS ===================
# cort_samples = ["BC28", "BC3", "BC9"]
# oil_samples = ["BC15", "BC14", "BC13"]
sample_key = 'sample'
groupby = 'subclass_name'
condition_key = 'group'
sc.settings.figdir = figs_dir

tf_regulons = data_dir / 'RAT_TF_regulons_df.csv' #installed previously on login node
net = pd.read_csv(tf_regulons)

# =================== INPUT ===================
#NOTE this object contains min 10 cells per sample so were good

# slide_tags_merged = glob.glob(f'{data_dir}/*.h5ad')
# adata_path = Path(slide_tags_merged[0])
# print(f'Using {adata_path} as the file path ...')
# adata = sc.read_h5ad(adata_path)
# adata.obs.loc[adata.obs[sample_key].isin(cort_samples), 'group'] = 'CORT'
# adata.obs.loc[adata.obs[sample_key].isin(oil_samples), 'group'] = 'OIL'

# ================ AFTER creating the DE object ================
# dea_path = output_base / 'first_try' / 'dea_result.csv'
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
lr_res.to_csv(lr_res_path)

# # ======================= Intracellular signaling networks =======================
# # utily function to select top n interactions
# def select_top_n(d, n=None):
#     d = dict(sorted(d.items(), key=lambda item: abs(item[1]), reverse=True))
#     return {k: v for i, (k, v) in enumerate(d.items()) if i < n}

# source_label = '017 CA3 Glut'
# target_label = '025 CA2-FC-IG Glut'

# # NOTE: We sort by the absolute value of the interaction stat
# lr_stats = lr_res[lr_res['source'].isin([source_label]) & lr_res['target'].isin([target_label])].copy()
# lr_stats = lr_stats.sort_values('interaction_stat', ascending=False, key=abs)

# # select receptors based on interaction stats
# lr_dict = lr_stats.set_index('receptor')['interaction_stat'].to_dict()
# input_scores = select_top_n(lr_dict, n=10)

# # First, let's transform the DEA statistics into a DF
# # we will use these to estimate deregulated TF activity
# dea_wide = dea_df[[groupby, 'stat']].reset_index(names='genes').pivot(index=groupby, columns='genes', values='stat')
# dea_wide = dea_wide.fillna(0)
# print(dea_wide)

# # Run Enrichment Analysis
# estimates, pvals = dc.mt.ulm(mat=dea_wide, net=net)
# estimates.T.sort_values(target_label, key=abs, ascending=False).head()

# # selecting the top TFs 
# tf_data = estimates.copy()
# tf_dict = tf_data.loc[target_label].to_dict()
# output_scores = select_top_n(tf_dict, n=5)

# # obtain ppi (protein-protein interaction) network
# ppis = op.interactions.OmniPath().get(genesymbols = True)

# ppis['mor'] = ppis['is_stimulation'].astype(int) - ppis['is_inhibition'].astype(int)
# ppis = ppis[(ppis['mor'] != 0) & (ppis['curation_effort'] >= 5) & ppis['consensus_direction']] 

# input_pkn = ppis[['source_genesymbol', 'mor', 'target_genesymbol']]
# input_pkn.columns = ['source', 'mor', 'target']
# input_pkn.head()

# # convert the PPI network into a knowledge graph
# prior_graph = li.mt.build_prior_network(input_pkn, input_scores, output_scores, verbose=True)

# # calculating node weights 
# temp = adata[adata.obs[groupby] == target_label].copy()

# node_weights = pd.DataFrame(temp.X.getnnz(axis=0) / temp.n_obs, index=temp.var_names)
# node_weights = node_weights.rename(columns={0: 'props'})
# node_weights = node_weights['props'].to_dict()

# # ===== CORNETO =====
# df_res, problem = li.mt.find_causalnet(
#     prior_graph, 
#     input_scores, 
#     output_scores, 
#     node_weights,
#     # penalize (max_penalty) nodes with counts in less than 0.1 of the cells
#     node_cutoff=0.1, 
#     max_penalty=1,
#     # the penaly of those in > 0.1 prop of cells set to:
#     min_penalty=0.01,
#     edge_penalty=0.1,
#     verbose=False,
#     max_runs=50, # NOTE that this repeats the solving either until the max runs are reached
#     stable_runs=10, # or until X number of consequitive stable runs are reached (i.e. no new edges are added)
#     solver='gurobi' # 'scipy' is available by default, but often results in suboptimal solutions
#     )

# cn.methods.carnival.visualize_network(df_res)