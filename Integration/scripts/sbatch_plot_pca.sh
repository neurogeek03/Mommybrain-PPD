#!/bin/bash
#SBATCH --job-name=plot_pca
#SBATCH --output=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --account=def-shreejoy
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

source /scratch/mfafouti/miniforge3/bin/activate
conda activate anndata_env
cd /scratch/mfafouti/Mommybrain/Integration/scripts
python plot_pca_interactive.py
