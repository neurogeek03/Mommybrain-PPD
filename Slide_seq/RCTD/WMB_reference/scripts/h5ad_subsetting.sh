#!/bin/bash
#SBATCH --job-name=split_h5ad  # Job name
#SBATCH --account=rrg-shreejoy          
#SBATCH --time=15:00:00                  # Wall time
#SBATCH --nodes=1
#SBATCH --ntasks=1                      # How many programs you are running
#SBATCH --output=/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/logs/%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/logs/%x_%j.err
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
conda activate anndata_env

cd /scratch/mfafouti/Mommybrain/Slide_seq/RCTD/WMB_reference/scripts

python split_h5ad.py

conda deactivate