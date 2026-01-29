#!/bin/bash

# ===== Usage =====
usage() {
    echo "Usage: $0 [--input-dir <input_directory>] [--output-dir <output_directory>]"
    exit 1
}

# ===== ARGUMENTS =====
ROOT="/scratch/mfafouti/Mommybrain/Slide_tags/Neuronchat"
INPUT_DIR="$ROOT/data"
OUTPUT_DIR=""

# Parse positional arguments                                                                                                          
if [ -z "$1" ]; then 
        usage                                                                                                                             
fi                                                                                                                                    
INPUT_DIR="$1"                                                                                                                        
OUTPUT_DIR="$2" # This will be empty if not provided  

cd $ROOT
echo "Using ROOT directory: $ROOT"

# ===== SAMPLES =====
SAMPLES=("BC13" "BC14" "BC28" "BC15" "BC3" "BC9")

# ===== DETERMINE OUTPUT DIRECTORY =====
if [ -z "$OUTPUT_DIR" ]; then
    read -p "Enter a base run name for this batch: " run_name
    base_folder_root="$ROOT/out/auto/$run_name"
else
    base_folder_root="$OUTPUT_DIR"
fi

echo "Base output directory will be: $base_folder_root"
mkdir -p "$base_folder_root"


# ===== MAIN LOOP =====
for SAMPLE in "${SAMPLES[@]}"; do
    echo "---------------------------------"
    echo "Processing sample: $SAMPLE"
    echo "---------------------------------"

    sample_base_folder="$base_folder_root/$SAMPLE"
    mkdir -p "$sample_base_folder"
    echo "Output for $SAMPLE will be saved to: $sample_base_folder"

    # ===== RUN PYTHON SCRIPT TO SUBSET =====
    source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
    conda activate sc_env

    echo "Subsetting sample: $SAMPLE"
    python3 "$ROOT/scripts/subset_sample.py" \
        --input "$INPUT_DIR" \
        --sample "$SAMPLE" \
        --output_base "$sample_base_folder"

    # ===== RUN R script through apptainer =====
    echo "Running NeuronChat for $SAMPLE"
    apptainer exec \
        --bind "$sample_base_folder":/mnt/data \
        --bind "$sample_base_folder":/mnt/out \
        --bind ./scripts:/mnt/scripts \
        /scratch/mfafouti/docker/neuronchat_full.sif \
        Rscript /mnt/scripts/test_neuronchat.R -d "/mnt/data" -s "$SAMPLE"
done

echo "All samples processed."
 