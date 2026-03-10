#!/bin/bash
#SBATCH --job-name=scanorama_test   # Job name
#SBATCH --account=rrg-shreejoy          
#SBATCH --time=5:00:00                  # Wall time
#SBATCH --nodes=1
#SBATCH --ntasks=1                      # How many programs you are running
#SBATCH --output=/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/classifier/logs/%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/classifier/logs/%x_%j.err
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# Activate conda environment
source /scratch/mfafouti/miniforge3/bin/activate

conda activate sc_env

# Go to project directory
cd /scratch/mfafouti/Mommybrain/Slide_seq/keons_single_cell/classifier/scripts

# Run your Python script
python label_transfer_scanorama.py