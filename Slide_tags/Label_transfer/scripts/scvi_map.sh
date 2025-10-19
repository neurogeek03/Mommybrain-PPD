#!/bin/bash
#SBATCH --job-name=Scvi_mapping         # Job name
#SBATCH --account=rrg-shreejoy  
#SBATCH --output=single_gpu_job_%j.out    # Output file (%j = job ID)
#SBATCH --nodes=1                         # Request 1 node
#SBATCH --gpus-per-node=1                 # Request 1 GPU
#SBATCH --output=/scratch/mfafouti/Mommybrain/Slide_tags/Label_transfer/logs/%x_%A_%a.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Slide_tags/Label_transfer/logs/%x_%A_%a.err
#SBATCH --time=00:10:00                  

source /scratch/mfafouti/miniforge3/etc/profile.d/conda.sh
conda activate sc_env

cd /scratch/mfafouti/Mommybrain/Slide_tags/Label_transfer/scripts

nvidia-smi

python scvi_label_transfer.py

