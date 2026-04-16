#!/bin/bash
# Run 04_plot_edgeR.py barplot --horizontal for all comparisons in an edger_lrt directory.
# Usage: bash scripts/run_barplots.sh <edger_lrt_dir>
# Example: bash scripts/run_barplots.sh /scratch/mfafouti/Mommybrain/Slide_tags/EdgeR/runs/seq_test/edger_lrt

set -euo pipefail

EDGER_DIR="${1:?Usage: bash run_barplots.sh <edger_lrt_dir>}"
SCRIPTS_DIR="$(dirname "$(realpath "$0")")"
SC_ENV="sc_env"

RUN_DIR="$(dirname "$(realpath "$EDGER_DIR")")"
OUT_DIR="$RUN_DIR/figures/SummaryBarPlots"
mkdir -p "$OUT_DIR"

for COMP_DIR in "$EDGER_DIR"/*/; do
  COMP=$(basename "$COMP_DIR")

  echo "Barplot: $COMP"
  conda run -n "$SC_ENV" python "$SCRIPTS_DIR/04_plot_edgeR.py" \
    --input-dir "$COMP_DIR" \
    barplot \
    --horizontal \
    --output "$OUT_DIR/${COMP}_barplot.png"
done

echo "Done. Figures saved under $OUT_DIR"
