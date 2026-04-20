#!/bin/bash
# Serve the static spatial gene expression viewer.
#
# Run export_binary.py first to generate viewer_data/.
# Usage: bash launch_static.sh
#
# Then on your LOCAL machine open a new terminal and run:
#   ssh -L 8080:<HOSTNAME>:8080 mfafouti@trillium.alliancecan.ca
# where <HOSTNAME> is printed below.
# Then open http://localhost:8080/viewer.html

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/viewer_data"
PORT=8080

if [ ! -f "$DATA_DIR/manifest.json" ]; then
  echo "ERROR: $DATA_DIR/manifest.json not found."
  echo "Run first: python export_binary.py --config test.conf"
  exit 1
fi

echo "============================================"
echo " Spatial Gene Expression Viewer (static)"
echo "============================================"
echo " Node:     $(hostname)"
echo " Serving:  $DATA_DIR"
echo " Port:     $PORT"
echo ""
echo " SSH tunnel command (run on your local machine):"
echo "   ssh -L ${PORT}:$(hostname):${PORT} mfafouti@trillium.alliancecan.ca"
echo ""
echo " Then open: http://localhost:${PORT}/viewer.html"
echo "============================================"

python -m http.server "$PORT" --directory "$DATA_DIR"
