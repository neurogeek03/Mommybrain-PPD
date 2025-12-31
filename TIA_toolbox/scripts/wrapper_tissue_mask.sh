#!/bin/bash

# ===== ARGUMENTS =====
# ROOT="$1"                     # optional: input dir root
RESOLUTION_TIA="$1"
SAMPLE_ARG="$2"               # optional: specific sample ID
ROOT="/scratch/mfafouti/Mommybrain/TIA_toolbox"

echo "Using ROOT: $ROOT"

# ===== FIXED INPUT PATH =====
INPUT_DIR="$ROOT/data/objects"
echo "Using INPUT_DIR: $INPUT_DIR"

# ===== PROMPT FOR RUN NAME =====
read -p "Enter run name explaining how you modified the pipeline this time: " run_name
base_folder="$ROOT/out/auto/$run_name"

# Define subfolders
masks="$base_folder/masks"
anndata_filtered="$base_folder/filtered_objects"
plots="$base_folder/spatial_plots"
tia_input="$base_folder/tia_input"

mkdir -p "$masks" "$anndata_filtered" "$plots" "$tia_input"

# ===== DETERMINE TARGET FILES =====
if [[ -n "$SAMPLE_ARG" ]]; then
    echo "Sample provided: $SAMPLE_ARG"
    FILES=("$INPUT_DIR"/"$SAMPLE_ARG"*.h5ad)
    if [[ ! -e "${FILES[0]}" ]]; then
        echo "No files found for sample '$SAMPLE_ARG' in $INPUT_DIR"
        exit 1
    fi
else
    echo "No sample provided → processing all files"
    FILES=("$INPUT_DIR"/*)
fi

# ===== MAIN LOOP =====
for f in "${FILES[@]}"; do
    echo "Processing: $f"
    sample=$(basename "$f" | cut -d'_' -f1)
    echo "$sample"

    source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
    conda activate anndata_env
    echo "Anndata environment activated!"

    # ==========================================================
    # 1. Plot UMI colored plots
    # ==========================================================
    echo "Plotting UMI colored plots ..."

    if ls "$tia_input"/"$sample"* 1>/dev/null 2>&1; then
        echo "→ Skipping UMI plots: sample '$sample' already exists in $tia_input"
    else
        python "$ROOT/scripts/plot_umi_brains.py" \
            -s "$sample" \
            -i "$f" \
            -o "$tia_input" \
            -p "$plots"
    fi

    # ==========================================================
    # 2. Run tissue segmentation & save plots + mask
    # ==========================================================
    module load apptainer/1.3.5
    APPTAINER_IMG="/scratch/mfafouti/docker/tiatoolbox_latest.sif"

    echo "Segmenting UMI-colored image for $sample ..."

    if ls "$masks"/"$sample"* 1>/dev/null 2>&1; then
        echo "→ Skipping segmentation: mask for '$sample' already exists in $masks"
    else
        apptainer exec "$APPTAINER_IMG" \
            python "$ROOT/scripts/tissue_mask_slideseq.py" \
                -s "$sample" \
                -i "$tia_input" \
                -o "$masks" \
                -p "$plots" \
                -res "$RESOLUTION_TIA"
    fi

    # ==========================================================
    # 3. Resize mask to fit the slide-seq puck & subset object
    # ==========================================================
    if ls "$anndata_filtered"/"$sample"* 1>/dev/null 2>&1; then
        echo "→ Skipping object subsetting: sample '$sample' already exists in $anndata_filtered"
    else
        python "$ROOT/scripts/resize_and_subset.py" \
            -s "$sample" \
            -i "$f" \
            -m "$masks" \
            -o "$anndata_filtered"
    fi

    echo "Finished sample $sample"
    echo "----------------------------------------"
    
done

    
