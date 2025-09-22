#!/bin/bash

MISSING_FILE="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/out_RCTD_all/missing_subsets.txt"
SBATCH_SCRIPT="/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/sbatch_RCTD_array.sh"

while IFS=: read -r SAMPLE SUBSETS; do
    if [ "$SUBSETS" != "None" ]; then
        for subset in $(echo $SUBSETS | tr ',' ' '); do
            echo "Submitting $SAMPLE subset $subset ..."
            sbatch --array=$subset \
                   --export=SAMPLE_ID=$SAMPLE,SEEKER_DIR=/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/new_RCTD_run/subsets_100/$SAMPLE \
                   $SBATCH_SCRIPT
        done
    fi
done < "$MISSING_FILE"