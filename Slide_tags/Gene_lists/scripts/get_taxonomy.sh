#!/bin/bash

# Base directory for ABC Atlas taxonomy metadata
BASE_DIR="/scratch/mfafouti/Mommybrain/Slide_tags/Gene_lists/out/data/abc_atlas_data/metadata/WMB-taxonomy/20230630"
mkdir -p "$BASE_DIR"
cd "$BASE_DIR"

# URL of the S3 folder
S3_URL="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/WMB-taxonomy/20230630/"

# Recursively download all files
wget -r -np -nH --cut-dirs=5 -R "index.html*" -c "$S3_URL"

echo "✅ All files downloaded recursively to $BASE_DIR"