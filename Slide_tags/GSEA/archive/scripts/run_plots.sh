#!/bin/bash
# ==========================================================================
# Title:        Plot gProfiler GO:BP results
# Description:  Runs 08_plot_gprofiler.py for all three cell type groups
# Author:       Maria Eleni Fafouti
# Usage:        bash scripts/run_plots.sh config/gsea_config.conf
# ==========================================================================
set -euo pipefail

GSEA_DIR="$(dirname "$(dirname "$(realpath "$0")")")"
SCRIPTS_DIR="$GSEA_DIR/scripts"

if [ $# -ge 1 ] && [ -f "$1" ]; then
  echo "Loading config: $1"
  source "$1"
else
  echo "ERROR: Please provide a config file."
  echo "Usage: bash scripts/run_plots.sh config/gsea_config.conf"
  exit 1
fi

GPROFILER_OUT="$GSEA_OUT_DIR/$RUN_NAME/gprofiler"

cd "$GSEA_DIR"
uv sync --quiet

echo "Plotting ${SOURCE:-GO:BP} results..."
uv run python "$SCRIPTS_DIR/03_plot_gprofiler.py" "$GPROFILER_OUT" "$GPROFILER_OUT/plots" \
  --top-n "${TOP_N:-30}" \
  --min-celltypes "${MIN_CELLTYPES:-2}" \
  --source "${SOURCE:-GO:BP}" \
  ${CELLTYPES_CSV:+--celltypes-csv "$CELLTYPES_CSV"}

echo "Done. Plots saved to: $GPROFILER_OUT/plots"
