#!/bin/bash
# Extract cell counts and sample metadata from h5ad.
# Runs extract_counts.py inside the sc_env conda environment.
#
# Usage:
#   bash scripts/extract_counts.sh [--h5ad PATH] [--out_dir PATH]
#
# Defaults match the main filtered object; override as needed.

set -euo pipefail

H5AD="/scratch/mfafouti/Mommybrain/Slide_tags/Filtering/out/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad"
OUT_DIR="out"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --h5ad)    H5AD="$2";    shift 2 ;;
        --out_dir) OUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

echo "H5AD:    $H5AD"
echo "OUT_DIR: $OUT_DIR"

conda run -n sc_env python scripts/extract_counts.py \
    --h5ad        "$H5AD" \
    --sample_col  "sample" \
    --celltype_col "subclass_name" \
    --treatment_col "treatment" \
    --out_dir     "$OUT_DIR"