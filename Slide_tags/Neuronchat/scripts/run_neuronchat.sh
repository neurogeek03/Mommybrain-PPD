#!/bin/bash

# ===== ARGUMENTS =====
ROOT="${1:-/scratch/mfafouti/Mommybrain/Slide_tags/Neuronchat}"

cd $ROOT
echo "Using ROOT directory: $ROOT"

# ===== PROMPT FOR RUN NAME =====
read -p "Enter run name explaining how you modified the pipeline this time: " run_name
base_folder="$ROOT/out/auto/$run_name"
mkdir -p "$base_folder"
echo "Output will be saved to: $base_folder"

# ===== RUN R script through apptainer =====
apptainer exec \
    --bind ./data:/mnt/data \
    --bind "$base_folder":/mnt/out \
    --bind ./scripts:/mnt/scripts \
    /scratch/mfafouti/docker/neuronchat_full.sif \
    Rscript /mnt/scripts/test_neuronchat.R 