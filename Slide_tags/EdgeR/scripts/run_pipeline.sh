#!/bin/bash
# ==========================================================================
# Title:        EdgeR DE Pipeline
# Description:  End-to-end pipeline: pseudobulk → EdgeR → plots → LIANA CSV
# Author:       Maria Eleni Fafouti
# Usage:        bash scripts/run_pipeline.sh  (from EdgeR/ directory)
# ==========================================================================
set -euo pipefail

# ========== CONFIGURATION ==========
WORKSPACE="/scratch/mfafouti/Mommybrain/Slide_tags/EdgeR"
SCRIPTS_DIR="$WORKSPACE/scripts"
APPTAINER_SIF="/scratch/mfafouti/docker/edger.sif"
SEURAT_ENV="/scratch/mfafouti/miniforge3/envs/seurat_env"
SC_ENV="sc_env"

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
fi

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

# ========== CREATE RUN FOLDER ==========
RUN_DIR="$WORKSPACE/runs/$RUN_NAME"

if [ -d "$RUN_DIR" ]; then
  echo "WARNING: '$RUN_DIR' already exists. Reusing existing folder."
fi

mkdir -p "$RUN_DIR/pseudobulk_outputs"
mkdir -p "$RUN_DIR/edger_lrt"
mkdir -p "$RUN_DIR/figures"
cp "$COMPARISONS_CSV" "$RUN_DIR/comparisons.csv"

echo "Run folder: $RUN_DIR"
echo ""

# ========== STEP 1: PSEUDOBULK ==========
echo "[1/5] Pseudobulking from AnnData..."
conda run -n "$SC_ENV" python "$SCRIPTS_DIR/01_pseudobulk_from_anndata.py" \
  --h5ad "$H5AD_FILE" \
  --output-dir "$RUN_DIR/pseudobulk_outputs" \
  --celltype-col "$CELLTYPE_COL" \
  --gene-symbol-col "$GENE_SYMBOL_COL"

# ========== STEP 2: EdgeR LRT (Apptainer) ==========
echo ""
echo "[2/5] Running EdgeR LRT (Apptainer)..."
apptainer exec \
  --bind "$WORKSPACE:/workspace" \
  --bind "$RUN_DIR/pseudobulk_outputs:/workspace/out/new_march_26/pseudobulk_outputs" \
  --bind "$RUN_DIR/edger_lrt:/workspace/out/edger_lrt" \
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
COMPARISONS=$(conda run -n "$SC_ENV" python "$_PARSE_SCRIPT" "$RUN_DIR/comparisons.csv")
rm "$_PARSE_SCRIPT"

for COMP in $COMPARISONS; do
  COMP_INPUT="$RUN_DIR/edger_lrt/$COMP"
  COMP_FIGURES="$RUN_DIR/figures/$COMP"
  mkdir -p "$COMP_FIGURES/volcanos"

  echo ""
  echo "--- Comparison: $COMP ---"

  echo "[3/5] Interactive volcano plots..."
  conda run -n "$SC_ENV" python "$SCRIPTS_DIR/03_plot_volcano.py" \
    --input-dir "$COMP_INPUT" \
    --output-dir "$COMP_FIGURES"

  echo "[4/5] Static volcano plots..."
  conda run -n "$SC_ENV" python "$SCRIPTS_DIR/04_plot_edgeR.py" \
    --input-dir "$COMP_INPUT" \
    volcanos --output-dir "$COMP_FIGURES/volcanos"

  echo "[5/5] LIANA DEA CSV..."
  conda run -n "$SC_ENV" python "$SCRIPTS_DIR/05_make_liana_dea_csv.py" \
    "$COMP_INPUT"

done

echo ""
echo "Pipeline complete. Results in: $RUN_DIR"
