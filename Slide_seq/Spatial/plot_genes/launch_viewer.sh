#!/bin/bash
# Launch the spatial gene expression Dash viewer.
#
# Usage: bash launch_viewer.sh
#
# Then on your LOCAL machine, open a new terminal and run:
#   ssh -L 8050:<HOSTNAME>:8050 mfafouti@trillium.alliancecan.ca
# where <HOSTNAME> is printed below when the script starts.
# Then open http://localhost:8050 in your browser.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/test.conf"
CONDA_ENV="sc_env"

echo "============================================"
echo " Spatial Gene Expression Viewer"
echo "============================================"
echo " Node:   $(hostname)"
echo " Config: $CONFIG"
echo ""
echo " SSH tunnel command (run on your local machine):"
echo "   ssh -L 8050:$(hostname):8050 mfafouti@trillium.alliancecan.ca"
echo ""
echo " Then open: http://localhost:8050"
echo "============================================"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

cd "$SCRIPT_DIR"
python main.py --config "$CONFIG"