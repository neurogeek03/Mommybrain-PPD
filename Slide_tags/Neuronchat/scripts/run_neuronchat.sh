#!/bin/bash

# ===== Usage =====
usage() {
    echo "Usage: $0 [--input-dir <input_directory>] [--output-dir <output_directory>] output dir is optional"
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
# SAMPLES=("BC13" "BC14" "BC28" "BC15" "BC3" "BC9")
SAMPLES=("CORT" "OIL")

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

    # ===== Step 1: Run Python Script to Subset =====
    echo "-----> Checking for: Subset CSVs"
    matrix_file="$sample_base_folder/${SAMPLE}_expression_matrix_cell_subclass.csv"
    meta_file="$sample_base_folder/${SAMPLE}_metadata.csv"

    if [ -f "$matrix_file" ] && [ -f "$meta_file" ]; then
        echo "Subset data for $SAMPLE found. Skipping."
    else
        echo "Subset data for $SAMPLE not found. Running subsetting script."
        source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
        conda activate sc_env

        python3 "$ROOT/scripts/subset_sample.py" \
            --input "$INPUT_DIR" \
            --sample "$SAMPLE" \
            --output_base "$sample_base_folder"
    fi
    echo "<----- Finished checking for: Subset CSVs"
    echo

    # ===== Step 2: Run NeuronChat R script =====
    echo "-----> Checking for: NeuronChat RDS object"
    nc_object_file="$sample_base_folder/${SAMPLE}_neuronchat_object.rds" # This is the object the plotting script needs

    if [ -f "$nc_object_file" ]; then
        echo "NeuronChat object for $SAMPLE found. Skipping analysis."
    else
        echo "NeuronChat object not found. Running analysis."
        apptainer exec \
            --bind "$sample_base_folder":/mnt/data \
            --bind "$sample_base_folder":/mnt/out \
            --bind ./scripts:/mnt/scripts \
            /scratch/mfafouti/docker/neuronchat_full.sif \
            Rscript /mnt/scripts/test_neuronchat.R --data_dir "/mnt/data" --sample "$SAMPLE"
    fi
    echo "<----- Finished checking for: NeuronChat RDS object"
    echo

    # ===== Step 3: Run Plotting Script =====
    echo "-----> Checking for: Output plots"
    plot_output_file="$sample_base_folder/figures/ranked_signaling_plots.pdf"

    if [ -f "$plot_output_file" ]; then
        echo "Plots for $SAMPLE found. Skipping plotting."
    else
        echo "Plots not found. Running plotting script."
        apptainer exec \
            --bind "$sample_base_folder":/mnt/data \
            --bind ./scripts:/mnt/scripts \
            /scratch/mfafouti/docker/neuronchat_full.sif \
            Rscript /mnt/scripts/plot_neuronchat.R --data_dir "/mnt/data" --sample "$SAMPLE"
    fi
    echo "<----- Finished checking for: Output plots"
    echo
done

echo "All samples processed."


# BASE_INPUT_DIR_HOST="out/auto/symbols_test"
# COMPARISON_OUTPUT_DIR_HOST="out/auto/symbols_test_comparison"


# nc_run_dir="symbols_test"

# mkdir -p "$COMPARISON_OUTPUT_DIR_HOST"


# apptainer shell \
#     --bind ./scripts:/mnt/scripts \
#     --bind ./out:/mnt/out \
#     /scratch/mfafouti/docker/neuronchat_full.sif 


# BASE_INPUT_DIR_CONTAINER="/mnt/out/auto/meta_added"
# COMPARISON_OUTPUT_DIR_CONTAINER="/mnt/out/auto/meta_test_comparison"

# Rscript /mnt/scripts/run_comparative_analysis.R --base_dir "$BASE_INPUT_DIR_CONTAINER" --output_dir "$COMPARISON_OUTPUT_DIR_CONTAINER"


# ########

# apptainer shell \
#     --bind ./out/auto/meta_added/BC:/mnt/input_object_dir \
#     --bind ./out/nmf_pathway_analysis/BC3_select_k:/mnt/output \
#     --bind ./scripts:/mnt/scripts \
#     /scratch/mfafouti/docker/neuronchat_full.sif 

#     Rscript /mnt/scripts/run_pathway_nmf_analysis.R --input_neuronchat_object_path "/mnt/input_object_dir/BC3_neuronchat_object.rds" --output_dir "/mnt/output" --k_range 2:10