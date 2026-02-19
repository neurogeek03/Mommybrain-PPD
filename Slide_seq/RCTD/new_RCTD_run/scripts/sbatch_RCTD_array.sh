#!/bin/bash
#SBATCH --job-name=RCTD_array_subsets_mouse   # Job name
#SBATCH --account=def-shreejoy          
#SBATCH --time=8:00:00                  # Wall time
#SBATCH --nodes=1
#SBATCH --ntasks=1                      # How many programs you are running
#SBATCH --mem=500G
#SBATCH --output=/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/logs/array_runs/%x_%A_%a.out
#SBATCH --error=/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/logs/array_runs/%x_%A_%a.err
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# ========== VALIDATE INPUT ==========
if [ -z "$SAMPLE_ID" ] || [ -z "$SEEKER_DIR" ]; then
  echo "❌ ERROR: SAMPLE_ID and SEEKER_DIR must be set via --export."
  exit 1
fi

# ========= SETUP LOG DIR ==========
LOG_BASE="/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/logs/array_runs"
LOG_DIR="${LOG_BASE}/${SAMPLE_ID}"
mkdir -p "$LOG_DIR"

# ========== MODULES ==========
module load apptainer

# ========== REFEERENCE ==========
# REF_RCTD="RATIFIED_ABC_ref_coronal_subclass_50k_25.rds" # CHANGE THIS ONLY
REF_RCTD="ABC_ref_coronal_subclass_50k_25.rds"

HOST_REF_PATH="/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/WMB_ref/${REF_RCTD}"
CONTAINER_REF_PATH="/workspace/WMB_ref/${REF_RCTD}"

# ========== VARIABLES ==========
# Derived container-visible path
CONTAINER_SEEKER_DIR="/workspace/subsets_100/${SAMPLE_ID}" #CHANGE FOLDER NAME FOR 2ND BATCH!!!!!!!!!!

# Get the matching file name from host
SEEKER_FILENAME=$(ls "$SEEKER_DIR"/${SAMPLE_ID}_subset_*.rds | sort -V | sed -n ${SLURM_ARRAY_TASK_ID}p)
SEEKER_FILENAME=$(basename "$SEEKER_FILENAME")

# Full container path to subsets
SEEKER_FILE="${CONTAINER_SEEKER_DIR}/${SEEKER_FILENAME}"

# Extract subset number from filename
SUBSET_NUM=$(echo "$SEEKER_FILENAME" | grep -oP "subset_\K[0-9]+")

# Create output dir
OUTPUT_DIR_HOST="/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/${SAMPLE_ID}"
mkdir -p "$OUTPUT_DIR_HOST"

OUTPUT_DIR="/workspace/out_RCTD_all/${SAMPLE_ID}"  # Container-visible path

# Output filenames
OUTPUT_FILENAME="${SAMPLE_ID}_subset_${SUBSET_NUM}_RCTD_results.csv"
OUTPUT_FILE="${OUTPUT_DIR}/${OUTPUT_FILENAME}"

OUTPUT_SEURAT_FILENAME="${SAMPLE_ID}_subset_${SUBSET_NUM}_seurat.rds"
OUTPUT_SEURAT="${OUTPUT_DIR}/${OUTPUT_SEURAT_FILENAME}"

OUTPUT_RCTD_FILENAME="${SAMPLE_ID}_subset_${SUBSET_NUM}_RCTD.rds"
OUTPUT_RCTD="${OUTPUT_DIR}/${OUTPUT_RCTD_FILENAME}"

OUTPUT_WEIGHTS_FILENAME="${SAMPLE_ID}_subset_${SUBSET_NUM}_weights.csv"
OUTPUT_WEIGHTS="${OUTPUT_DIR}/${OUTPUT_WEIGHTS_FILENAME}"

# ========== EXECUTE ==========
apptainer exec \
    --bind /scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run:/workspace \
    /scratch/mfafouti/seurat_rctd_lib_fix.sif \
    Rscript /workspace/scripts/RCTD_analysis_no_parall_for_array.R "$SEEKER_FILE" "$CONTAINER_REF_PATH" "$OUTPUT_FILE" "$OUTPUT_SEURAT" "$OUTPUT_RCTD" "$OUTPUT_WEIGHTS" 