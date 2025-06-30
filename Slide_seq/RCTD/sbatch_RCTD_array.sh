#!/bin/bash
#SBATCH --job-name=RCTD_array_subsets
#SBATCH --account=def-shreejoy
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --array=1-100
#SBATCH --output=/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/logs/%x_%j.out
#SBATCH --error=/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/logs/%x_%j.err
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# ========== MODULES ==========
module load apptainer

# ========== VALIDATE INPUT ==========
if [ -z "$SAMPLE_ID" ] || [ -z "$SEEKER_DIR" ]; then
  echo "❌ ERROR: SAMPLE_ID and SEEKER_DIR must be set via --export."
  exit 1
fi

# ========== VARIABLES ==========
# Derived container-visible path
CONTAINER_SEEKER_DIR="/workspace/subsets_100_batch1/${SAMPLE_ID}" #CHANGE FOLDER NAME FOR 2ND BATCH!!!!!!!!!!
CONTAINER_REF_PATH="/workspace/No_bad_celltype_ref_slide_tags.rds"

# Get the matching file name from host
SEEKER_FILENAME=$(ls "$SEEKER_DIR"/${SAMPLE_ID}_subset_*.rds | sort -V | sed -n ${SLURM_ARRAY_TASK_ID}p)
SEEKER_FILENAME=$(basename "$SEEKER_FILENAME")

# Full container path to subset
SEEKER_FILE="${CONTAINER_SEEKER_DIR}/${SEEKER_FILENAME}"

# Extract subset number from filename
SUBSET_NUM=$(echo "$SEEKER_FILENAME" | grep -oP "subset_\K[0-9]+")

# Create output dir
OUTPUT_DIR_HOST="/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/out_RCTD_all/${SAMPLE_ID}"
mkdir -p "$OUTPUT_DIR_HOST"

OUTPUT_DIR="/workspace/out_RCTD_all/${SAMPLE_ID}"  # Container-visible path

# Output filenames
OUTPUT_FILENAME="${SAMPLE_ID}_subset_${SUBSET_NUM}_RCTD_results.csv"
OUTPUT_FILE="${OUTPUT_DIR}/${OUTPUT_FILENAME}"

OUTPUT_SEURAT_FILENAME="${SAMPLE_ID}_subset_${SUBSET_NUM}_seurat.rds"
OUTPUT_SEURAT="${OUTPUT_DIR}/${OUTPUT_SEURAT_FILENAME}"

OUTPUT_RCTD_FILENAME="${SAMPLE_ID}_subset_${SUBSET_NUM}_RCTD.rds"
OUTPUT_RCTD="${OUTPUT_DIR}/${OUTPUT_RCTD_FILENAME}"

# ========== EXECUTE ==========
apptainer exec \
  --bind /home/mfafouti/miniforge3/envs/rctd_env:/opt/rctd_env \
  --bind /home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run:/workspace \
  --bind /home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/subsets_100_batch1:/workspace/subsets_100_batch1 \
  --bind /home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/out_RCTD_all:/workspace/out_RCTD_all \
  --bind /home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/No_bad_celltype_ref_slide_tags.rds:/workspace/No_bad_celltype_ref_slide_tags.rds \
  /home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/seurat_v5.sif \
  Rscript /workspace/RCTD_analysis_no_parall_for_array.R "$SEEKER_FILE" "$CONTAINER_REF_PATH" "$OUTPUT_FILE" "$OUTPUT_SEURAT" "$OUTPUT_RCTD"