import anndata as ad
from pathlib import Path
from cell_type_mapper.cli.from_specified_markers import main as from_specified_markers_main
from cell_type_mapper.taxonomy.taxonomy_tree import TaxonomyTree
from cell_type_mapper.cli.precompute_stats_abc import (PrecomputationABCRunner)
import json

scratch_dir = Path("/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/out/data/scratch")

# Defining dirs 
precompute_dir = scratch_dir / 'precompute'
reference_dir = scratch_dir / 'reference'
query_dir = scratch_dir / 'query'

# Point to your local copy of WMB taxonomy metadata
taxonomy_dir = Path("/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/out/data/abc_atlas_data/metadata/WMB-taxonomy")

# Adjust extensions to match what you actually have (csv/parquet)
cluster_annotation_path = taxonomy_dir / "cluster_annotation_term.csv"
cluster_membership_path = taxonomy_dir / "cluster_to_cluster_annotation_membership.csv"

precompute_config = {
    'output_path': str(precompute_dir / 'precomputed_stats.h5'),
    'hierarchy': ['CCN20230722_CLAS',
                  'CCN20230722_SUBC',
                  'CCN20230722_SUPT',
                  'CCN20230722_CLUS'],
    'h5ad_path_list': h5ad_path_list,
    'cell_metadata_path': str(training_set_path),
    'cluster_annotation_path': str(cluster_annotation_path),
    'cluster_membership_path': str(cluster_membership_path),
    'n_processors': 4,
    'split_by_dataset': True,
    'do_pruning': True,
    'tmp_dir': str(scratch_dir),
    'clobber': True
}

precomputation_runner = PrecomputationABCRunner(
    args=[],
    input_data=precompute_config
)

precomputation_runner.run()

# # === SETUP ===
# # Use project folder for writable directory
# scratch_dir = Path("/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/scratch")
# scratch_dir.mkdir(parents=True, exist_ok=True)

# # Your 5 h5ad input files
# input_files = [
#     "BC13_collapsed_mouse_genes.h5ad",
#     "BC14_collapsed_mouse_genes.h5ad",
#     "BC15_collapsed_mouse_genes.h5ad",
#     "BC28_collapsed_mouse_genes.h5ad",
#     "BC3_collapsed_mouse_genes.h5ad",
#     "BC9_collapsed_mouse_genes.h5ad",
# ]
# # Option 1: map separately

# # Paths to Allen Brain Atlas reference files (downloaded earlier)
# marker_path = scratch_dir / "mouse_markers_230821.json"
# precompute_path = scratch_dir / "precomputed_stats_ABC_revision_230821.h5"

# # === DOWNLOAD AND PRINT TAXONOMY TREE ===
# precomputed_path = scratch_dir / "precomputed_stats_ABC_revision_230821.h5"
# tree = TaxonomyTree.from_precomputed_stats_h5(precomputed_path)
# print("Levels:", tree.hierarchy)

# # === DEFINE REGIONS TO DROP ===
# # You need to know their "taxonomy level" and node name (class or region ID)
# # For example, to exclude cerebellum, olfactory areas, and hypothalamus:
# regions_to_exclude = [
#     ('region', 'Cerebellum'),
#     ('region', 'Hypothalamus'),
#     ('region', 'Olfactory areas')
# ]

# # Or if you’re excluding taxonomy classes instead:
# # regions_to_exclude = [('class', 'CS20230722_CLAS_17'), ...]

# # === LOOP OVER FILES ===
# for h5ad_path in input_files:
#     output_json = scratch_dir / f"{Path(h5ad_path).stem}_mapping.json"
#     output_csv = scratch_dir / f"{Path(h5ad_path).stem}_mapping.csv"

#     mapping_config = {
#         'query_path': str(h5ad_path),
#         'extended_result_path': str(output_json),
#         'csv_result_path': str(output_csv),
#         'tmp_dir': str(scratch_dir),
#         'max_gb': 10,
#         'cloud_safe': False,
#         'verbose_stdout': False,
#         'type_assignment': {
#             'normalization': 'raw',
#             'n_processors': 4,
#             'chunk_size': 10000,
#             'bootstrap_iteration': 100,
#             'bootstrap_factor': 0.5,
#             'rng_seed': 42
#         },
#         'precomputed_stats': {
#             'path': str(precompute_path)
#         },
#         'query_markers': {
#             'serialized_lookup': str(marker_path)
#         },
#         'nodes_to_drop': regions_to_exclude,
#         'drop_level': 'CCN20230722_SUPT'  # taxonomy level name — check your taxonomy file
#     }

#     print(f"=== Running mapping for {h5ad_path} ===")
#     runner = FromSpecifiedMarkersRunner(args=[], input_data=mapping_config)
#     runner.run()
#     print(f"Finished {h5ad_path}")