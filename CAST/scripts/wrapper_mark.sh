#!/usr/bin/env bash
set -euo pipefail

# Args
# Get the sample ID from the first argument to the bash script
SAMPLE_ID="${1:-B01}"  # defaults to B01 if no argument is provided

# ================= CAST MARK - INPUTS =================
# Activate conda env
source /scratch/miniforge3/etc/profile.d/conda.sh
conda activate anndata_env

# Step 1: run Python script
python scripts/cast_mark_input.py

# Step 2: Apptainer step
SCRATCH="/scratch/mfafouti"
PROJECT="$SCRATCH/Mommybrain/CAST"
APPTAINER_IMG="$SCRATCH/docker/cast_first.sif"
OUTPUT_DIR="$PROJECT/out"

cd PROJECT

apptainer exec "$APPTAINER_IMG" \
    python cast_mark_run.py \
    --input "$INPUT_DIR" \
    --output "$OUTPUT_DIR"
