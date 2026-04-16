#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../config/pipeline.conf"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found at $CONFIG"
    exit 1
fi

source "$CONFIG"

# Prompt for run name and create output directory
read -rp "Enter a name for this run: " RUN_NAME
RUN_DIR="${SCRIPT_DIR}/../out/runs/${RUN_NAME}"

START_FROM=0
if [[ -d "$RUN_DIR" ]]; then
    echo "Run directory already exists: $RUN_DIR"
    echo "  0) Start fresh (re-run all steps)"
    echo "  2) Resume from Step 02  (skip preprocessing)"
    echo "  3) Resume from Step 03  (plots only)"
    read -rp "Start from step [0/2/3]: " START_FROM
fi

mkdir -p "$RUN_DIR"
echo "Output will be written to: $RUN_DIR"

echo "============================================================"
echo " LIANA CCC Pipeline"
echo " Config: $CONFIG"
echo " Run: $RUN_NAME  (starting from step ${START_FROM})"
echo "============================================================"

PYTHON="conda run -n liana python"
cd "$SCRIPT_DIR"

# -------- Step 00: Preprocess --------
if [[ $START_FROM -le 0 ]]; then
    echo ""
    echo "--- Step 00: Preprocess AnnData ---"

    MERGE_ARGS=()
    if [[ ${#MERGE_CELLTYPES[@]} -gt 0 ]]; then
        MERGE_ARGS+=(--merge_celltypes "${MERGE_CELLTYPES[@]}")
        MERGE_ARGS+=(--merge_label "$MERGE_LABEL")
    fi

    $PYTHON 00_preprocess_adata.py \
        --cort_samples    "${CORT_SAMPLES[@]}" \
        --oil_samples     "${OIL_SAMPLES[@]}" \
        --sample_key      "$SAMPLE_KEY" \
        --celltype_key    "$CELLTYPE_KEY" \
        --condition_key   "$CONDITION_KEY" \
        --n_top_genes     "$N_TOP_GENES" \
        --run_dir         "$RUN_DIR" \
        "${MERGE_ARGS[@]}"
else
    echo ""
    echo "--- Step 00: Preprocess AnnData  [skipped] ---"
fi

# -------- Step 02: LIANA df_to_lr --------
if [[ $START_FROM -le 2 ]]; then
    echo ""
    echo "--- Step 02: LIANA df_to_lr (EdgeR) ---"

    $PYTHON 02_liana_w_edgeR.py \
        --dea_path      "$DEA_PATH" \
        --celltype_key  "$CELLTYPE_KEY" \
        --condition_key "$CONDITION_KEY" \
        --resource_name "$RESOURCE_NAME" \
        --expr_prop     "$EXPR_PROP" \
        --run_dir       "$RUN_DIR"
else
    echo ""
    echo "--- Step 02: LIANA df_to_lr (EdgeR)  [skipped] ---"
fi

# -------- Step 03: Plots --------
if [[ $START_FROM -le 3 ]]; then
    echo ""
    echo "--- Step 03: Plot cell type pairs ---"

    $PYTHON 03_plot_celltype_pairs.py \
        --p_val   "$P_VAL" \
        --run_dir "$RUN_DIR"
else
    echo ""
    echo "--- Step 03: Plot cell type pairs  [skipped] ---"
fi

# -------- Step 04: TF enrichment --------
if [[ $START_FROM -le 4 ]]; then
    echo ""
    echo "--- Step 04: TF enrichment (dc.mt.ulm) ---"

    $PYTHON 04_intracellular_signaling.py \
        --dea_path     "$DEA_PATH" \
        --celltype_key "$CELLTYPE_KEY" \
        --run_dir      "$RUN_DIR"
else
    echo ""
    echo "--- Step 04: TF enrichment  [skipped] ---"
fi

echo ""
echo "============================================================"
echo " Pipeline complete."
echo "============================================================"
