#!/bin/bash
#SBATCH --job-name=00_prepare_data
#SBATCH --output=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH --time=1:00:00
#SBATCH --account=def-shreejoy
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

source /scratch/mfafouti/miniforge3/bin/activate
conda activate anndata_env
cd /scratch/mfafouti/Mommybrain/Integration/scripts
python 00_prepare_data.py
