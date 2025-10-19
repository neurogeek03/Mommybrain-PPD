#!/bin/bash

# Base directory for ABC Atlas metadata
BASE_DIR="/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/out/data/abc_atlas_data/metadata/WMB-10X/20230630"
mkdir -p "$BASE_DIR"
cd "$BASE_DIR"

# List of metadata files (add all files you need here)
metadata_files=(
    "cluster_annotations.csv"
    "cell_types.csv"
    "feature_annotations.csv"
    "sample_annotations.csv"
    "gene_metadata.csv"
    # Add any other files you see in the S3 folder
)

# Loop over each file and download
for file in "${metadata_files[@]}"; do
    url="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/WMB-10X/20230630/${file}"
    echo "⬇️ Downloading $file"
    wget -c "$url" -O "$file"
done

echo "✅ All metadata files downloaded to $BASE_DIR"