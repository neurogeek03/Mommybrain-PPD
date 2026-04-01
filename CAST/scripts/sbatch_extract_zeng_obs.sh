#!/bin/bash
#SBATCH --job-name=extract_zeng_obs
#SBATCH --output=/scratch/mfafouti/Mommybrain/CAST/reference/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/CAST/reference/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=6:00:00
#SBATCH --account=def-shreejoy
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# Activate conda environment
source /scratch/mfafouti/miniforge3/bin/activate

conda activate anndata_env

# Go to project directory
cd /scratch/mfafouti/Mommybrain/CAST/scripts

# Run Python script
python extract_zeng_obs.py