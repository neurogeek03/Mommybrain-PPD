#!/bin/bash
# ==========================================================================
# Title:        EdgeR DE Pipeline
# Description:  End-to-end pipeline: pseudobulk → EdgeR → plots → LIANA CSV
# Author:       Maria Eleni Fafouti
# Usage:        bash scripts/run_pipeline.sh  (from EdgeR/ directory)
# ==========================================================================
set -euo pipefail

# ========== CONFIGURATION ==========
WORKSPACE="/scratch/mfafouti/Mommybrain-PPD/Slide_tags/EdgeR"
SCRIPTS_DIR="$WORKSPACE/scripts"
APPTAINER_SIF="/scratch/mfafouti/docker/edger.sif"
SEURAT_ENV="/scratch/mfafouti/miniforge3/envs/seurat_env"
PYTHON="$WORKSPACE/.venv/bin/python"

# ========== PROMPT FOR INPUTS ==========
echo "=== EdgeR Pipeline Setup ==="
echo ""

if [ $# -ge 1 ] && [ -f "$1" ]; then
  echo "Loading config: $1"
  source "$1"
else
  read -rp "Run name: " RUN_NAME
  read -rp "Path to .h5ad file: " H5AD_FILE
  read -rp "Path to comparisons CSV: " COMPARISONS_CSV
  read -rp "Cell type column name [subclass_name]: " CELLTYPE_COL
  CELLTYPE_COL="${CELLTYPE_COL:-subclass_name}"
  read -rp "Gene symbol column name [gene_symbols]: " GENE_SYMBOL_COL
  GENE_SYMBOL_COL="${GENE_SYMBOL_COL:-gene_symbols}"
  read -rp "EdgeR model [lrt/qlf, default: lrt]: " EDGER_MODEL
  EDGER_MODEL="${EDGER_MODEL:-lrt}"
fi

EDGER_MODEL="${EDGER_MODEL:-lrt}"
LIB_SIZE_MIN="${LIB_SIZE_MIN:-50000}"
FILTER_MIN_COUNT="${FILTER_MIN_COUNT:-10}"
FILTER_MIN_TOTAL_COUNT="${FILTER_MIN_TOTAL_COUNT:-15}"
FILTER_LARGE_N="${FILTER_LARGE_N:-10}"
FILTER_MIN_PROP="${FILTER_MIN_PROP:-0.7}"
EDGER_SUBDIR="edger_${EDGER_MODEL}"

H5AD_FILE=$(realpath "$H5AD_FILE")
if [ ! -f "$H5AD_FILE" ]; then
  echo "ERROR: h5ad file not found: $H5AD_FILE"
  exit 1
fi

COMPARISONS_CSV=$(realpath "$COMPARISONS_CSV")
if [ ! -f "$COMPARISONS_CSV" ]; then
  echo "ERROR: comparisons CSV not found: $COMPARISONS_CSV"
  exit 1
fi

echo ""

# ========== CREATE OUTPUT FOLDERS ==========
PSEUDOBULK_DIR="$WORKSPACE/out/pseudobulk_outputs"
OUTPUT_BASE_DIR="${OUTPUT_BASE_DIR:-$WORKSPACE/out/runs}"
RUN_DIR="$OUTPUT_BASE_DIR/$RUN_NAME"

if [ -d "$RUN_DIR" ]; then
  echo "WARNING: '$RUN_DIR' already exists. Reusing existing folder."
fi

mkdir -p "$PSEUDOBULK_DIR"
mkdir -p "$RUN_DIR/$EDGER_SUBDIR"
mkdir -p "$RUN_DIR/figures"
cp "$COMPARISONS_CSV" "$RUN_DIR/comparisons.csv"

# Copy config into run folder for reproducibility
if [ $# -ge 1 ] && [ -f "$1" ]; then
  cp "$1" "$RUN_DIR/run_config.conf"
fi

# ========== LOGGING ==========
LOG_FILE="$RUN_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "Log file: $LOG_FILE"

echo "Pseudobulk dir: $PSEUDOBULK_DIR"
echo "Run folder:     $RUN_DIR"
echo ""

# ========== STEP 1: PSEUDOBULK ==========
EXISTING_COUNTS=$(find "$PSEUDOBULK_DIR" -name "*_counts.tsv" 2>/dev/null | wc -l)
if [ "$EXISTING_COUNTS" -gt 0 ]; then
  echo "[1/5] Skipping pseudobulk — $EXISTING_COUNTS count files already exist in $PSEUDOBULK_DIR"
else
  echo "[1/5] Pseudobulking from AnnData..."
  "$PYTHON" "$SCRIPTS_DIR/01_pseudobulk_from_anndata.py" \
    --h5ad "$H5AD_FILE" \
    --output-dir "$PSEUDOBULK_DIR" \
    --celltype-col "$CELLTYPE_COL" \
    --gene-symbol-col "$GENE_SYMBOL_COL"
fi

# ========== STEP 2: EdgeR LRT (Apptainer) ==========
echo ""
echo "[2/5] Running EdgeR (model: $EDGER_MODEL) (Apptainer)..."
apptainer exec \
  --env EDGER_MODEL="$EDGER_MODEL" \
  --env LIB_SIZE_MIN="$LIB_SIZE_MIN" \
  --env FILTER_MIN_COUNT="$FILTER_MIN_COUNT" \
  --env FILTER_MIN_TOTAL_COUNT="$FILTER_MIN_TOTAL_COUNT" \
  --env FILTER_LARGE_N="$FILTER_LARGE_N" \
  --env FILTER_MIN_PROP="$FILTER_MIN_PROP" \
  --bind "$WORKSPACE:/workspace" \
  --bind "$PSEUDOBULK_DIR:/workspace/out/pseudobulk_outputs" \
  --bind "$RUN_DIR/$EDGER_SUBDIR:/workspace/out/edger_out" \
  --bind "$RUN_DIR/comparisons.csv:/workspace/comparisons.csv" \
  --bind "$WORKSPACE/rostral_caudal.csv:/workspace/rostral_caudal.csv" \
  --bind "$SEURAT_ENV:/opt/seurat_env" \
  "$APPTAINER_SIF" \
  Rscript /workspace/scripts/02_run_edger.R

# ========== STEPS 3–5: Per comparison ==========
# Parse GroupA_vs_GroupB names from comparisons.csv
_PARSE_SCRIPT=$(mktemp /tmp/parse_comps_XXXXXX.py)
cat > "$_PARSE_SCRIPT" << 'PYEOF'
import pandas as pd, sys
df = pd.read_csv(sys.argv[1])
for _, row in df.iterrows():
    print(f"{row['GroupA_name']}_vs_{row['GroupB_name']}")
PYEOF
COMPARISONS=$("$PYTHON" "$_PARSE_SCRIPT" "$RUN_DIR/comparisons.csv")
rm "$_PARSE_SCRIPT"

for COMP in $COMPARISONS; do
  COMP_INPUT="$RUN_DIR/$EDGER_SUBDIR/$COMP"
  COMP_FIGURES="$RUN_DIR/figures/$COMP"
  mkdir -p "$COMP_FIGURES/volcanos"

  echo ""
  echo "--- Comparison: $COMP ---"

  echo "[3/5] Interactive volcano plots..."
  "$PYTHON" "$SCRIPTS_DIR/03_plot_volcano.py" \
    --input-dir "$COMP_INPUT" \
    --output-dir "$COMP_FIGURES"

  echo "[4/5] Static volcano plots..."
  "$PYTHON" "$SCRIPTS_DIR/04_plot_edgeR.py" \
    --input-dir "$COMP_INPUT" \
    volcanos --output-dir "$COMP_FIGURES/volcanos"

  echo "[5/5] LIANA DEA CSV..."
  "$PYTHON" "$SCRIPTS_DIR/05_make_liana_dea_csv.py" \
    "$COMP_INPUT"

done

# ========== BARPLOTS: All comparisons ==========
echo ""
echo "[+] Summary barplots..."
BARPLOT_OUT="$RUN_DIR/figures/SummaryBarPlots"
mkdir -p "$BARPLOT_OUT"
for COMP_DIR in "$RUN_DIR/$EDGER_SUBDIR"/*/; do
  COMP=$(basename "$COMP_DIR")
  echo "  Barplot: $COMP"
  "$PYTHON" "$SCRIPTS_DIR/04_plot_edgeR.py" \
    --input-dir "$COMP_DIR" \
    barplot \
    --horizontal \
    --output "$BARPLOT_OUT/${COMP}_barplot.png"
done

echo ""
echo "Pipeline complete. Results in: $RUN_DIR"
