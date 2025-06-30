#!/bin/bash

# Top-level directory containing one folder per sample (e.g., B08, B09, ...)
BASE_DIR="/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/subsets_100_batch1"

# Loop through each subdirectory
for sample_dir in "$BASE_DIR"/*; do
    # Only process directories
    if [ -d "$sample_dir" ]; then
        # Extract the sample name from the folder path (e.g., B08)
        sample_name=$(basename "$sample_dir")

        # Count how many subset_*.rds files exist in this sample folder
        num_subsets=$(ls "$sample_dir"/*_subset_*.rds 2>/dev/null | wc -l)

        # Skip if no subsets found
        if [ "$num_subsets" -eq 0 ]; then
            echo "⚠️  No subset files found for $sample_name, skipping..."
            continue
        fi

        echo "🟢 Submitting job for $sample_name with $num_subsets subsets..."

        # Submit the array job
        sbatch --array=1-${num_subsets} \
               --export=SAMPLE_ID="$sample_name",SEEKER_DIR="$sample_dir" \
               sbatch_RCTD_array.sh
    fi
done