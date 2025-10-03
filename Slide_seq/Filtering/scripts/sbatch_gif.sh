#!/bin/bash
#SBATCH --job-name=subset_RCTD_query   # Job name
#SBATCH --account=rrg-shreejoy          
#SBATCH --time=2:00:00                  # Wall time
#SBATCH --nodes=1
#SBATCH --ntasks=1                      # How many programs you are running
#SBATCH --output=/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/logs/%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/new_RCTD_run/logs/%x_%j.err
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
conda activate anndata_env 
echo "Anndata environment activated!"

cd /scratch/mfafouti/Mommybrain-PPD/Slide_seq/RCTD/Filtering/scripts

python gif_umap.py

conda deactivate 