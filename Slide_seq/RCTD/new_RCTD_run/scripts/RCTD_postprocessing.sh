#!/bin/bash
# ===== ARGUMENTS =====
INPUT_DIR="$1"   # first argument: input directory for merge_rctd_out.py
DEFAULT="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all"

# Use provided input or default
INPUT_DIR="${INPUT_DIR:-$DEFAULT}"
echo "Using INPUT_DIR: $INPUT_DIR"

# ===== FIXED BASE PATH =====
ROOT="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run"

# ===== PROMPT FOR RUN NAME =====
read -p "Enter run name explaining how you modified the pipeline this time: " run_name
base_folder="$ROOT/OUT/$run_name"

# Define subfolders
merged="$base_folder/merged_metadata_csvs"
anndata="$base_folder/anndata_objects"
plots="$base_folder/spatial_plots"
umaps="$base_folder/umaps"

mkdir -p "$merged" "$anndata" "$plots" "$umaps"

# Activate correct conda env 
source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
conda activate anndata_env 
echo "Anndata environment activated!"

# Run scripts in sequence, only if needed
echo "Merging metadata files from subsets..."
if [ -z "$(ls -A "$merged")" ]; then
    python $ROOT/scripts/merge_rctd_out.py -i "$INPUT_DIR" -o "$merged"
else
    echo "→ Skipping: merged metadata already exists in $merged"
fi

echo "Adding merged metadata to AnnData objects..."
if [ -z "$(ls -A "$anndata")" ]; then
    python $ROOT/scripts/add_rctd_meta_mouse.py -i "$merged" -o "$anndata"
else
    echo "→ Skipping: AnnData objects already exist in $anndata"
fi

echo "Plotting spatial data..."
if [ -z "$(ls -A "$plots")" ]; then
    python $ROOT/scripts/indiv_spatial_plots.py -i "$anndata" -o "$plots"
else
    echo "→ Skipping: spatial plots already exist in $plots"
fi

echo "Plotting UMAPs..."
if [ -z "$(ls -A "$umaps")" ]; then
    python $ROOT/scripts/run_umap.py -i "$anndata" -o "$umaps"
else
    echo "→ Skipping: UMAPs already exist in $umaps"
fi

echo "All scripts finished! Results in $base_folder/"

