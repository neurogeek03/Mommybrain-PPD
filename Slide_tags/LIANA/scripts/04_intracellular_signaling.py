import argparse
from datetime import datetime
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import decoupler as dc

parser = argparse.ArgumentParser(description='TF enrichment from EdgeR DEA using dc.mt.ulm')
parser.add_argument('--dea_path',     type=str, required=True,
                    help='Path to EdgeR DEA CSV (same as used in step 02)')
parser.add_argument('--celltype_key', type=str, default='subclass_name')
parser.add_argument('--run_dir',      type=str, required=True,
                    help='Path to the run output directory (e.g. out/runs/my_run)')
parser.add_argument('--tmin',         type=int, default=5,
                    help='Min number of regulon targets per TF (dc.mt.ulm tmin, default 5)')
parser.add_argument('--n_top_tfs',    type=int, default=15,
                    help='Number of top TFs per cell type to include in summary plot')
args = parser.parse_args()

stamp = datetime.now().strftime("%M%S%H_%Y%m%d")
print(f"------ Script started at {stamp} ------")

# Paths
project_path = Path.cwd().parents[0]
run_dir  = Path(args.run_dir)
out_dir  = run_dir / 'intracellular'
out_dir.mkdir(exist_ok=True, parents=True)
figs_dir = out_dir / 'figures'
figs_dir.mkdir(exist_ok=True, parents=True)

# =================== INPUT ===================

# DEA stats from EdgeR
dea_df = pd.read_csv(args.dea_path, index_col=0)
print(f"DEA loaded: {dea_df.shape[0]} gene-celltype rows, "
      f"{dea_df[args.celltype_key].nunique()} cell types")

# TF regulons — keep only the three columns dc.mt.ulm needs
tf_regulons_path = project_path / 'data' / 'RAT_TF_regulons_df.csv'
net = pd.read_csv(tf_regulons_path)[['source', 'target', 'weight']]
print(f"TF regulons loaded: {net['source'].nunique()} TFs, "
      f"{net['target'].nunique()} target genes, {len(net)} edges")

# lr_res from step 02
lr_res_path = run_dir / 'liana_edgeR' / 'lr_res.csv'
lr_res = pd.read_csv(lr_res_path, index_col=0)
print(f"lr_res loaded: {lr_res.shape[0]} LR pairs, "
      f"{lr_res['target'].nunique()} target cell types")

# =================== BUILD dea_wide ===================
# Reshape to (cell_type x gene) matrix of EdgeR stat values for dc.mt.ulm.
# Genes not tested in a cell type are filled with 0 (neutral, no enrichment signal).

dea_wide = (
    dea_df[[args.celltype_key, 'stat']]
    .reset_index(names='genes')
    .pivot(index=args.celltype_key, columns='genes', values='stat')
    .fillna(0)
)
print(f"dea_wide shape: {dea_wide.shape}  (cell types x genes)")

# =================== TF ENRICHMENT ===================
print(f"\nRunning dc.mt.ulm (tmin={args.tmin}) ...")

estimates, pvals = dc.mt.ulm(dea_wide, net, tmin=args.tmin)
print(f"TF activity estimates: {estimates.shape}  (cell types x TFs)")

# Save
estimates.to_csv(out_dir / 'tf_estimates.csv')
pvals.to_csv(out_dir / 'tf_pvals.csv')
print(f"Saved estimates -> {out_dir / 'tf_estimates.csv'}")
print(f"Saved p-values  -> {out_dir / 'tf_pvals.csv'}")

# =================== SUMMARY: top TFs per cell type ===================
records = []
for celltype in estimates.index:
    top = estimates.loc[celltype].abs().nlargest(args.n_top_tfs)
    for tf, _ in top.items():
        records.append({
            'cell_type': celltype,
            'TF': tf,
            'activity': estimates.loc[celltype, tf],
            'pval': pvals.loc[celltype, tf],
        })

top_tfs_df = pd.DataFrame(records)
top_tfs_df.to_csv(out_dir / 'top_tfs_per_celltype.csv', index=False)
print(f"Top {args.n_top_tfs} TFs per cell type saved.")

# =================== PLOT: TF activity heatmap ===================
# Show TFs that appear in the top-n list of at least one cell type.
top_tf_names = top_tfs_df['TF'].unique()
heatmap_mat  = estimates[top_tf_names]

fig, ax = plt.subplots(figsize=(min(1.0 * len(top_tf_names), 40), 0.4 * len(estimates)))
sns.heatmap(
    heatmap_mat,
    cmap='RdBu_r',
    center=0,
    linewidths=0.3,
    ax=ax,
    cbar_kws={'label': 'TF activity (ULM score)'},
)
ax.set_title('TF activity across cell types (dc.mt.ulm on EdgeR stats)')
ax.set_xlabel('Transcription Factor')
ax.set_ylabel('Cell type')
plt.tight_layout()
heatmap_path = figs_dir / 'tf_activity_heatmap.png'
fig.savefig(heatmap_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Heatmap saved -> {heatmap_path}")

# =================== PLOT: per-cell-type TF activity barplots ===================
# One barplot per cell type, saved to figures/per_celltype/.
per_ct_dir = figs_dir / 'per_celltype'
per_ct_dir.mkdir(exist_ok=True)

for celltype in estimates.index:
    scores  = estimates.loc[celltype].dropna()
    top_idx = scores.abs().nlargest(args.n_top_tfs).index
    plot    = scores[top_idx].sort_values()

    fig, ax = plt.subplots(figsize=(6, 0.4 * len(plot)))
    colors = ['#d62728' if v > 0 else '#1f77b4' for v in plot]
    ax.barh(plot.index, plot.values, color=colors)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('TF activity (ULM score)')
    ax.set_title(f'Top {args.n_top_tfs} TFs -- {celltype}')
    plt.tight_layout()

    # sanitise cell type name for use as filename
    safe_name = celltype.replace('/', '-').replace(' ', '_')
    fig.savefig(per_ct_dir / f'tf_activity_{safe_name}.png', dpi=150, bbox_inches='tight')
    plt.close()

print(f"Per-cell-type barplots saved -> {per_ct_dir}")

print(f"\n------ Script complete at {datetime.now().strftime('%M%S%H_%Y%m%d')} ------")