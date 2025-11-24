#!/bin/bash
#SBATCH --job-name=Cast_mark_input
#SBATCH --output=/scratch/mfafouti/Mommybrain/CAST/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/CAST/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --time=3:00:00
#SBATCH --account=def-shreejoy_cpu
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# Activate conda environment
source /scratch/mfafouti/miniforge3/bin/activate

conda activate anndata_env

# Go to project directory
cd /scratch/mfafouti/Mommybrain/CAST/scripts

# Run your Python script
python cast_mark_input.py