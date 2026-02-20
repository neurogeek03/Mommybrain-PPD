#!/bin/bash

# Usage: ./run_extraction_pipeline.sh <input_rds_path> <output_dir>
# Example: ./run_extraction_pipeline.sh out_RCTD_all/B01/B01_subset_1_RCTD.rds out_RCTD_all/csvs_from_rctd_obj

# Default values if no arguments provided
INPUT_RDS=${1:-"/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/B01/B01_subset_1_RCTD.rds"}
OUT_DIR=${2:-"/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/csvs_from_rctd_obj"}

# Paths
CONTAINER="/scratch/mfafouti/docker/seurat_rctd_lib_fix.sif"
R_SCRIPT="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/scripts/export_rctd_weights.R"
PY_SCRIPT="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/scripts/visualize_rctd_weights.py"

# Ensure output directory exists
mkdir -p "$OUT_DIR"

echo "--------------------------------------------------"
echo "STEP 1: Extracting weights from RDS inside Apptainer"
echo "Input: $INPUT_RDS"
echo "--------------------------------------------------"

apptainer exec "$CONTAINER" Rscript "$R_SCRIPT" "$INPUT_RDS" "$OUT_DIR"

echo ""
echo "--------------------------------------------------"
echo "STEP 2: Generating Heatmap Visualization"
echo "--------------------------------------------------"

# Assuming you are in the 'base' conda environment as seen in your terminal
python "$PY_SCRIPT" 
    "$OUT_DIR/rctd_full_weights.csv" 
    "$OUT_DIR/rctd_weights_heatmap.png"

echo ""
echo "Process Complete."
echo "CSV files and Heatmap are located in: $OUT_DIR"
echo "--------------------------------------------------"
