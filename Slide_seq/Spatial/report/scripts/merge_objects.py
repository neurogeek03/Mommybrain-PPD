import anndata as ad

neurons_path = '/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/out/semi_filtered_neurons/181663_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad'
non_neurons_path = '/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/out/broader_class_filter_NN/singlet_score_330_260k_cells/267456_umap_filtered_0_NEW_genelist_slide_seq_15.h5ad'

celltype_col = 'RCTD_first_type_rat'

out_path = '/scratch/mfafouti/Mommybrain/Slide_seq/Spatial/report/out_NN_seq/merged_neurons_non_neurons_slide_seq_15.h5ad'

neurons = ad.read_h5ad(neurons_path)
non_neurons = ad.read_h5ad(non_neurons_path)

merged = ad.concat(
    [neurons, non_neurons],
    join='outer',
    merge='same',
    label='source',
    keys=['neurons', 'non_neurons']
)

merged.write_h5ad(out_path)
print(f"Saved merged object ({merged.n_obs} cells x {merged.n_vars} genes) to:\n{out_path}")