#!/bin/bash
#SBATCH --job-name=preview
#SBATCH --output=/scratch/mfafouti/Mommybrain/Integration/tmp/preview_%j.out
#SBATCH --error=/scratch/mfafouti/Mommybrain/Integration/tmp/preview_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --account=def-shreejoy

source /scratch/mfafouti/miniforge3/bin/activate
conda activate anndata_env
python ~/preview/preview_anndata.py "/scratch/mfafouti/Mommybrain/Integration/out/configB/integrated.h5ad"
