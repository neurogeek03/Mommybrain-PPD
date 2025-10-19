#!/bin/bash
#SBATCH --job-name=Extract_anndata_files   # Job name
#SBATCH --account=rrg-shreejoy          
#SBATCH --time=2:00:00                  # Wall time
#SBATCH --nodes=1
#SBATCH --ntasks=1                      # How many programs you are running
#SBATCH --output=/scratch/mfafouti/Mommybrain/Slide_tags/Label_transfer/logs/%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Slide_tags/Label_transfer/logs/%x_%j.err
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
conda activate sc_env

cd /scratch/mfafouti/Mommybrain/Slide_tags/Label_transfer/scripts

python extract_anndata.py

conda deactivate