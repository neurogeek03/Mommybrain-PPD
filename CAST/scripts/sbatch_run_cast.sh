#!/bin/bash
#SBATCH --job-name=Cast_run_CPU
#SBATCH --output=/scratch/mfafouti/Mommybrain/CAST/logs/_%x_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/CAST/logs/_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --time=1:00:00
#SBATCH --account=def-shreejoy_cpu
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# Go to project directory
cd /scratch/mfafouti/Mommybrain/CAST/scripts

module load apptainer/1.3.5
# module load cuda/12.2

# echo "=== GPU visible to this job ==="
# nvidia-smi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "=== Running CAST inside container ==="
apptainer exec /scratch/mfafouti/docker/cast_first.sif python run_cast.py