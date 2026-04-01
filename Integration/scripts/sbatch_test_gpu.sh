#!/bin/bash
#SBATCH --job-name=test_gpu
#SBATCH --output=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Integration/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_2g.20gb:1
#SBATCH --account=def-shreejoy
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=END,FAIL

cd /scratch/mfafouti/Mommybrain/Integration
uv run python scripts/test_gpu.py
