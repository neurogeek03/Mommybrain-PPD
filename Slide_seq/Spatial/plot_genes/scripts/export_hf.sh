#!/bin/bash
#SBATCH --job-name=spatial_hf_export
#SBATCH --account=rrg-shreejoy
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/mfafouti/Mommybrain-PPD/Slide_seq/Spatial/plot_genes/logs/%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain-PPD/Slide_seq/Spatial/plot_genes/logs/%x_%j.err
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
conda activate sc_env

cd /scratch/mfafouti/Mommybrain-PPD/Slide_seq/Spatial/plot_genes

echo "[$(date)] Starting full-transcriptome binary export"

python export_binary.py --config config/hf_export.conf

echo "[$(date)] Export complete"

conda deactivate
