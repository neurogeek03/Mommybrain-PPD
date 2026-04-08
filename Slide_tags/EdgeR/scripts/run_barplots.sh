#!/bin/bash
# Run 04_plot_edgeR.py (barplot or volcanos) for all comparisons in an edger_lrt directory.
# Usage: bash scripts/run_barplots.sh <edger_lrt_dir> barplot|volcanos [conda_env]
# Example: bash scripts/run_barplots.sh /scratch/mfafouti/Mommybrain/Slide_tags/EdgeR/runs/seq_test/edger_lrt barplot sc_env

set -euo pipefail

EDGER_DIR="${1:?Usage: bash run_barplots.sh <edger_lrt_dir> barplot|volcanos [conda_env]}"
PLOT_TYPE="${2:?Usage: bash run_barplots.sh <edger_lrt_dir> barplot|volcanos [conda_env]}"
SC_ENV="${3:-sc_env}"
SCRIPTS_DIR="$(dirname "$(realpath "$0")")"

RUN_DIR="$(dirname "$(realpath "$EDGER_DIR")")"

for COMP_DIR in "$EDGER_DIR"/*/; do
  COMP=$(basename "$COMP_DIR")

  if [[ "$PLOT_TYPE" == "barplot" ]]; then
    OUT_DIR="$RUN_DIR/figures/SummaryBarPlots"
    mkdir -p "$OUT_DIR"
    echo "Barplot: $COMP"
    conda run -n "$SC_ENV" python "$SCRIPTS_DIR/04_plot_edgeR.py" \
      --input-dir "$COMP_DIR" \
      barplot \
      --horizontal \
      --output "$OUT_DIR/${COMP}_barplot.png" \
      --export-celltypes out/celltype_list.csv

  elif [[ "$PLOT_TYPE" == "volcanos" ]]; then
    OUT_DIR="$RUN_DIR/figures/Volcanos/$COMP"
    mkdir -p "$OUT_DIR"
    echo "Volcanos: $COMP"
    conda run -n "$SC_ENV" python "$SCRIPTS_DIR/04_plot_edgeR.py" \
      --input-dir "$COMP_DIR" \
      volcanos \
      --output-dir "$OUT_DIR"

  else
    echo "Unknown plot type: $PLOT_TYPE. Must be 'barplot' or 'volcanos'."
    exit 1
  fi
done

echo "Done. Figures saved under $RUN_DIR/figures/"