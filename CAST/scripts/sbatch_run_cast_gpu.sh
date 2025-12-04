#!/bin/bash
#SBATCH --job-name=Cast_STACK_test
#SBATCH --output=/scratch/mfafouti/Mommybrain/CAST/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/CAST/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=h100:1
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --account=def-shreejoy_cpu
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# Go to project directory
cd /scratch/mfafouti/Mommybrain/CAST/scripts

module load apptainer/1.3.5
# module load cuda/12.2

echo "=== GPU visible to this job ==="
nvidia-smi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:128

echo "=== Running CAST inside container ==="
apptainer exec --nv /scratch/mfafouti/docker/cast_first.sif python3 -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))" python 3_cast_stack.py