#!/bin/bash

# This script runs the pathway NMF analysis using NeuronChat's R script within an Apptainer container.
# It takes the full path to a NeuronChat object (.rds file) as an argument,
# creates a corresponding output directory, and executes the analysis.

# Check if an argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <path_to_neuronchat_object.rds>"
  echo "Example: $0 /path/to/my_neuronchat_object.rds"
  exit 1
fi

NEURONCHAT_OBJECT_PATH="$1"

# Validate if the provided path is an .rds file
if [[ ! "$NEURONCHAT_OBJECT_PATH" =~ \.rds$ ]]; then
  echo "Error: The provided path must be an .rds file."
  exit 1
fi

# Ensure the NeuronChat object file exists
if [ ! -f "$NEURONCHAT_OBJECT_PATH" ]; then
  echo "Error: NeuronChat object file not found at '$NEURONCHAT_OBJECT_PATH'"
  exit 1
fi

# Get the directory containing the NeuronChat object
INPUT_OBJECT_DIR_HOST=$(dirname "${NEURONCHAT_OBJECT_PATH}")
INPUT_OBJECT_FILENAME=$(basename "${NEURONCHAT_OBJECT_PATH}")

# Define the base output directory relative to the current working directory
BASE_OUTPUT_DIR="./out/nmf_pathway_analysis"

# Extract a unique identifier from the filename for the output directory
# For example, if INPUT_OBJECT_FILENAME is "BC3_neuronchat_object.rds", ID will be "BC3"
# If the filename does not contain an underscore, it will use the full filename without extension
ID=$(echo "${INPUT_OBJECT_FILENAME}" | sed -E 's/(_neuronchat_object)?\.rds$//')

OUTPUT_DIR_HOST="${BASE_OUTPUT_DIR}/${ID}_select_k"

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR_HOST}"

# Define the path to the scripts directory, relative to where this script is located
# Assuming this script is in ./scripts/
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# Path to the Apptainer image
APPTAINER_IMAGE="/scratch/mfafouti/docker/neuronchat_full.sif" # Adjust if the sif file is located elsewhere

echo "Starting NMF Pathway Analysis for: ${NEURONCHAT_OBJECT_PATH}"
echo "Output will be saved to: ${OUTPUT_DIR_HOST}"
echo "Using Apptainer image: ${APPTAINER_IMAGE}"

# Execute the Apptainer command
# We bind the host directories to specific mount points within the container
# /mnt/input_object_dir will contain the input .rds file
# /mnt/output will be the directory where results are written
# /mnt/scripts will contain the R scripts
apptainer exec \
  --bind "${INPUT_OBJECT_DIR_HOST}:/mnt/input_object_dir" \
  --bind "${OUTPUT_DIR_HOST}:/mnt/output" \
  --bind "${SCRIPT_DIR}:/mnt/scripts" \
  "${APPTAINER_IMAGE}" \
  Rscript /mnt/scripts/run_pathway_nmf_analysis.R \
  --input_neuronchat_object_path "/mnt/input_object_dir/${INPUT_OBJECT_FILENAME}" \
  --output_dir "/mnt/output" \
  --k 3

if [ $? -eq 0 ]; then
  echo "NMF Pathway Analysis completed successfully."
else
  echo "NMF Pathway Analysis failed."
  exit 1
fi
