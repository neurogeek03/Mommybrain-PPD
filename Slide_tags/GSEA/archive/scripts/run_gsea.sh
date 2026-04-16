#!/bin/bash
# ==========================================================================
# Title:        GSEA Pipeline
# Description:  Ranked gene lists → gProfiler enrichment per cell type
# Author:       Maria Eleni Fafouti
# Usage:        bash scripts/run_gsea.sh config/gsea_config.conf
# ==========================================================================
set -euo pipefail

GSEA_DIR="$(dirname "$(dirname "$(realpath "$0")")")"
SCRIPTS_DIR="$GSEA_DIR/scripts"

# ========== LOAD CONFIG ==========
if [ $# -ge 1 ] && [ -f "$1" ]; then
  echo "Loading config: $1"
  source "$1"
else
  echo "ERROR: Please provide a config file."
  echo "Usage: bash scripts/run_gsea.sh config/gsea_config.conf"
  exit 1
fi

RNK_OUT="$GSEA_OUT_DIR/$RUN_NAME/rnk"
GPROFILER_OUT="$GSEA_OUT_DIR/$RUN_NAME/gprofiler"

# ========== SETUP ENV ==========
echo "Setting up uv environment..."
cd "$GSEA_DIR"
uv sync
echo ""

# ========== STEP 1: RANKED GENE LISTS ==========
if [ -d "$RNK_OUT" ] && [ -n "$(ls "$RNK_OUT"/*.rnk 2>/dev/null)" ]; then
  echo "[1/2] Skipping ranking — .rnk files already exist in $RNK_OUT"
else
  echo "[1/2] Making ranked gene lists (.rnk)..."
  uv run python "$SCRIPTS_DIR/01_make_rnk.py" "$EDGER_INPUT_DIR" "$RNK_OUT"
fi

# ========== STEP 2: GPROFILER ==========
echo ""
echo "[2/2] Running gProfiler..."
uv run python "$SCRIPTS_DIR/02_run_gprofiler.py" "$RNK_OUT" "$GPROFILER_OUT" \
  --organism "$ORGANISM" \
  --mode "${GPROFILER_MODE:-ordered}" \
  --edger-dir "$EDGER_INPUT_DIR" \
  --de-fdr "${DE_FDR:-0.05}" \
  --background-mode "${BACKGROUND_MODE:-global}" \
  --gprofiler-fdr "${GPROFILER_FDR:-0.05}"

echo ""
echo "Done. Results in: $GSEA_OUT_DIR/$RUN_NAME"
echo "To plot: bash scripts/run_plots.sh config/gsea_config.conf"
