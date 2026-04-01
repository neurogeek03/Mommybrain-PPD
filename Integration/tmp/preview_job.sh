#!/bin/bash
#SBATCH --job-name=preview
#SBATCH --output=/scratch/mfafouti/Mommybrain/Integration/tmp/preview_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Integration/tmp/preview_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --account=def-shreejoy

source /scratch/mfafouti/miniforge3/bin/activate
conda activate anndata_env

echo "========== Slide-tags (n=207903) =========="
python ~/preview/preview_anndata.py "/scratch/mfafouti/Mommybrain/Integration/data/Merged_n=207903_all_metadata_slide_tags.h5ad"

echo ""
echo "========== Slide-seq neurons =========="
python ~/preview/preview_anndata.py "/scratch/mfafouti/Mommybrain/Integration/data/adata_neuron_subset_181663_singlets_score_330.h5ad"

echo ""
echo "========== Slide-seq non-neurons =========="
python ~/preview/preview_anndata.py "/scratch/mfafouti/Mommybrain/Integration/data/adata_nn_subset_78689_singlets_score_330.h5ad"
