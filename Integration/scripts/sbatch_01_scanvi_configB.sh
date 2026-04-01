#!/bin/bash
#SBATCH --job-name=01_scanvi_configB
#SBATCH --output=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_2g.20gb:2
#SBATCH --account=def-shreejoy
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

cd /scratch/mfafouti/Mommybrain/Integration

# Full run (train scVI + scANVI):
#   uv run python scripts/01_scanvi.py config/01_scanvi_configB.conf
# Resume from saved scVI model (skip scVI, run scANVI only):
uv run python scripts/01_scanvi.py config/01_scanvi_configB.conf --resume
