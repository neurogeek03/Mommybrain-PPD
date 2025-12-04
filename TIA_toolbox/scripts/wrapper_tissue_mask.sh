#!/bin/bash

# ===== ARGUMENTS =====
ROOT="$1"   # first argument: input directory
DEFAULT="/scratch/mfafouti/Mommybrain/TIA_toolbox"

# Use provided input or default
ROOT="${ROOT:-$DEFAULT}"     # FIX: ROOT_DIR → ROOT
echo "Using ROOT: $ROOT"

# ===== FIXED INPUT PATH =====
INPUT_DIR="$ROOT/data/objects"      # FIX: use correct path
echo "Using INPUT_DIR: $INPUT_DIR"

# ===== PROMPT FOR RUN NAME =====
read -p "Enter run name explaining how you modified the pipeline this time: " run_name
base_folder="$ROOT/out/auto/$run_name"

# Define subfolders
masks="$base_folder/masks"
anndata_filtered="$base_folder/filtered_objects"
plots="$base_folder/spatial_plots"
tia_input="$base_folder/tia_input"

mkdir -p "$masks" "$anndata_filtered" "$plots" "$tia_input"

# ===== Looping over all adata objects =====
# FIX: remove duplicated path; just loop inside INPUT_DIR
for f in "$INPUT_DIR"/*; do
    echo "Processing: $f"
    sample=$(basename "$file" | cut -d'_' -f1)
    echo "$sample"

    source /scratch/miniforge3/etc/profile.d/conda.sh
    conda activate anndata_env
    echo "Anndata environment activated!"

    # ===== 1. Plot UMI colored plots =====
    echo "Plotting UMI colored plots ..."

    # FIX: test for empty merged directory; merge only if empty
    if [ -z "$(ls -A "$tia_input" 2>/dev/null)" ]; then
        python "$ROOT/scripts/plot_umi_brains.py" -s "$sample" -i "$f" -o "$tia_input" -p "$plots"
    else
        echo "→ Skipping: merged metadata already exists in $tia_input"
    fi

done

    # ===== 2. Run tissue segmentation & save diagnostic plots + mask
    # APPTAINER_IMG="/scratch/mfafouti/docker/tiatoolbox_latest.sif"

    # apptainer exec "$APPTAINER_IMG" \
    #     python cast_mark_run.py \
    #     --input "$tia_input" \
    #     --output "$OUTPUT_DIR"