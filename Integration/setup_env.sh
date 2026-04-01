#!/bin/bash
# Set up integration environment with uv
# Run on a compute node (needs internet access or pre-cached packages)

ENV_NAME="integration_env"
ENV_PATH="/scratch/mfafouti/envs/${ENV_NAME}"

# Create venv with uv
uv venv "${ENV_PATH}" --python 3.11

# Install packages
uv pip install --python "${ENV_PATH}/bin/python" \
    "scvi-tools>=1.2" \
    "scanpy>=1.10" \
    "scib>=1.1" \
    "plotly>=5.0" \
    "matplotlib" \
    "pyyaml" \
    "leidenalg" \
    "igraph"

# PyTorch with CUDA (scvi-tools may pull CPU-only by default)
# Uncomment if GPU torch is not installed automatically:
# uv pip install --python "${ENV_PATH}/bin/python" \
#     torch --index-url https://download.pytorch.org/whl/cu121

echo ""
echo "Environment created at: ${ENV_PATH}"
echo "Activate with: source ${ENV_PATH}/bin/activate"
