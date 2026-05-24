#!/bin/bash
# Run this after the SLURM export job completes.
# Requires: huggingface_hub installed in base conda env
# Usage: bash scripts/upload_to_hf.sh

set -euo pipefail

HF_USER="marlenfafouti"
DATASET_REPO="${HF_USER}/mommybrain-ppd-spatial"
SPACE_REPO="${HF_USER}/spatial-gene-viewer"
DATA_DIR="/scratch/mfafouti/Mommybrain-PPD/Slide_seq/Spatial/plot_genes/viewer_data_hf"
VIEWER_HTML="/scratch/mfafouti/Mommybrain-PPD/Slide_seq/Spatial/plot_genes/viewer_hf.html"

echo "=== Step 1: Login to Hugging Face ==="
echo "You will be prompted for your HF token."
echo "Get it from: https://huggingface.co/settings/tokens (read+write)"
huggingface-cli login

echo ""
echo "=== Step 2: Create dataset repo (skip if already exists) ==="
huggingface-cli repo create "${DATASET_REPO}" --type dataset --private || echo "Repo may already exist, continuing..."

echo ""
echo "=== Step 3a: Upload manifest + sample bins ==="
huggingface-cli upload "${DATASET_REPO}" "${DATA_DIR}/manifest.json" manifest.json \
    --repo-type dataset \
    --commit-message "Add manifest and sample bins"
huggingface-cli upload "${DATASET_REPO}" "${DATA_DIR}/B03.bin" B03.bin --repo-type dataset
huggingface-cli upload "${DATASET_REPO}" "${DATA_DIR}/B14.bin" B14.bin --repo-type dataset
huggingface-cli upload "${DATASET_REPO}" "${DATA_DIR}/B41.bin" B41.bin --repo-type dataset

echo ""
echo "=== Step 3b: Upload gene bins in batches (avoids rate limit) ==="
echo "This will take a while (~12 GB across many small commits)..."
huggingface-cli upload "${DATASET_REPO}" "${DATA_DIR}/genes" genes \
    --repo-type dataset \
    --every 500 \
    --num-workers 1 \
    --commit-message "Upload gene expression bins"

echo ""
echo "=== Step 4: Create Space repo (skip if already exists) ==="
huggingface-cli repo create "${SPACE_REPO}" --type space --space-sdk static --private || echo "Space may already exist, continuing..."

echo ""
echo "=== Step 5: Upload viewer HTML to Space ==="
huggingface-cli upload "${SPACE_REPO}" "${VIEWER_HTML}" index.html \
    --repo-type space \
    --commit-message "Add spatial gene expression viewer"

echo ""
echo "=== Done ==="
echo "Dataset : https://huggingface.co/datasets/${DATASET_REPO}"
echo "Viewer  : https://huggingface.co/spaces/${SPACE_REPO}"
