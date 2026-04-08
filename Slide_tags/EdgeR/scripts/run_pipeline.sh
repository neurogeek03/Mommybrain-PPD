#!/bin/bash
# ==========================================================================
# Title:        EdgeR DE Pipeline
# Description:  End-to-end pipeline: pseudobulk → EdgeR → plots → LIANA CSV
# Author:       Maria Eleni Fafouti
# Usage:        bash scripts/run_pipeline.sh  (from EdgeR/ directory)
# ==========================================================================
set -euo pipefail

# ========== LOAD CONFIG ==========
echo "=== DE analysis using EdgeR LRT ==="
echo "... loading config"
echo ""

if [ $# -ge 1 ] && [ -f "$1" ]; then
  echo "Loading config: $1"
  source "$1"
else
  echo "ERROR: Please provide a config file."
  echo "Usage: bash scripts/run_pipeline.sh config/run_config.conf"
  exit 1
fi

SCRIPTS_DIR="$WORKSPACE/scripts"

MERGE_CELLTYPES=("038 DG-PIR Ex IMN" "045 OB-STR-CTX Inh IMN")
MERGE_LABEL="038-045 IMN"

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
_PSEUDOBULK_ARGS=(
  --h5ad "$H5AD_FILE"
  --output-dir "$RUN_DIR/pseudobulk_outputs"
  --celltype-col "$CELLTYPE_COL"
  --gene-symbol-col "$GENE_SYMBOL_COL"
)
_PSEUDOBULK_ARGS+=(--merge-celltypes "${MERGE_CELLTYPES[@]}" --merge-label "$MERGE_LABEL")
conda run -n "$SC_ENV" python "$SCRIPTS_DIR/01_pseudobulk_from_anndata.py" \
  "${_PSEUDOBULK_ARGS[@]}"

# ========== STEP 2: EdgeR LRT (Apptainer) ==========
echo ""
echo "[2/5] Running EdgeR LRT (Apptainer)..."
source /cvmfs/soft.computecanada.ca/custom/software/lmod/lmod/init/bash
module load apptainer
apptainer exec \
  --bind "$WORKSPACE:/workspace" \
  --bind "$RUN_DIR/pseudobulk_outputs:/workspace/out/new_march_26/pseudobulk_outputs" \
  --bind "$RUN_DIR/edger_lrt:/workspace/out/edger_lrt" \
  --bind "$RUN_DIR/comparisons.csv:/workspace/comparisons.csv" \
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

  echo "[5/6] LIANA DEA CSV..."
  conda run -n "$SC_ENV" python "$SCRIPTS_DIR/05_make_liana_dea_csv.py" \
    "$COMP_INPUT"

  echo "[6/6] GSEA ranked gene lists (.rnk)..."
  conda run -n "$SC_ENV" python "$SCRIPTS_DIR/06_make_rnk.py" \
    --input-dir "$COMP_INPUT" \
    --output-dir "$RUN_DIR/gsea/$COMP"

done

echo ""
echo "Pipeline complete. Results in: $RUN_DIR"


#  conda run -n anndata_env python scripts/06_make_rnk.py --input-dir runs/merge_imn/edger_lrt/OIL_vs_CORT --output-dir runs/merge_imn/gsea/OIL_vs_CORT