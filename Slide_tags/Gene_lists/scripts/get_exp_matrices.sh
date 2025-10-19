#!/bin/bash

# Base directory for ABC Atlas data
BASE_DIR="/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/out/data/abc_atlas_data/expression_matrices/WMB-10Xv3/20230630"
mkdir -p "$BASE_DIR"
cd "$BASE_DIR"

# List of feature matrices
feature_matrices=(
    "WMB-10Xv3-HPF" "WMB-10Xv3-Isocortex-1" "WMB-10Xv3-PAL"
    "WMB-10Xv3-STR" "WMB-10Xv3-CTXsp" "WMB-10Xv3-HY"
    "WMB-10Xv3-OLF" "WMB-10Xv3-TH" "WMB-10Xv3-MB" "WMB-10Xv3-Isocortex-2"
)

# Loop over each feature matrix and download the raw H5AD file
for feature in "${feature_matrices[@]}"; do
    url="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/${feature}-raw.h5ad"
    echo "⬇️ Downloading $feature"
    wget -c "$url" -O "${feature}-raw.h5ad"
done

echo "✅ All files downloaded to $BASE_DIR"

METADATA_DIR="/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/out/data/abc_atlas_data/metadata/WMB-10X/20230630"
mkdir -p "$METADATA_DIR"
cd "$METADATA_DIR"

# List of metadata files
metadata_files=(
    "gene.csv"
    "region_of_interest_metadata.csv"
    "cell_metadata.csv"
)

for file in "${metadata_files[@]}"; do
    url="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/WMB-10X/20230630/${file}"
    echo "⬇️ Downloading $feature"
    wget -c "$url" -O "${feature}-raw.h5ad"
done