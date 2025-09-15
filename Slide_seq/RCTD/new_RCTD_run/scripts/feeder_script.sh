#!/bin/bash
#SBATCH --job-name=feeder
#SBATCH --time=12:00:00
#SBATCH --mem=2G

# Base directory where your sample folders live
BASE_DIR="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/subsets_100"

# Explicit list of samples
samples=("B11" "B14" "B18" "B19" "B21" "B31" "B41" "B42")

for sample in "${samples[@]}"; do
    # Wait until you have fewer than 50 jobs in the queue
    while [ $(squeue -u $USER | wc -l) -ge 300 ]; do
        echo "[$(date)] Queue full (≥50 jobs). Sleeping for 2 hours..."
        sleep 7200   # check every 2 hours
    done

    echo "[$(date)] Submitting job for sample: $sample"

    cd /scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run

    # Submit the job array for this sample
    sbatch --job-name=RCTD_${sample} \
           --array=1-200 \
           --export=SAMPLE_ID=$sample,SEEKER_DIR=$BASE_DIR/$sample \
           sbatch_RCTD_array.sh
done
