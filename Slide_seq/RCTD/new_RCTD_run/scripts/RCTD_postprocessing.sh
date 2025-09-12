#!/bin/bash

# ===== FIXED BASE PATH =====
ROOT="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/results"

# ===== PROMPT FOR RUN NAME =====
read -p "Enter run name explaining how you modified the pipeline this time: " run_name
base_folder="$ROOT/$run_name"

# Define subfolders
merged="$ROOT/$base_folder/merged_metadata_csvs"
anndata="$ROOT/$base_folder/anndata_objects"
plots="$ROOT/$base_folder/spatial_plots"

mkdir -p "$merged" "$anndata" "$plots"

# Run scripts in sequence
echo "Merging metadata files from subsets..."
python merge_rctd_out.py -o "$merged"

echo "Adding merged metadata to AnnData objects..."
python add_rctd_meta_mouse.py -i "$merged" -o "$anndata"

echo "Plotting spatial data..."
python indiv_spatial_plots.py -i "$anndata" -o "$plots"

echo "All scripts finished! Results in $base_folder/"

