#!/bin/bash
#SBATCH --job-name=Seurat_Integration   # Name of your job
#SBATCH --account=def-shreejoy                
#SBATCH --time=3:00:00                        
#SBATCH --nodes=1                           
#SBATCH --ntasks=1                             
#SBATCH --cpus-per-task=8                  
#SBATCH --mem=128G                              
#SBATCH --output=/home/mfafouti/scratch/Mommybrain_marlen/Slide_tags/Single_cell_analysis/Integration/Integration_out/%x_%j.out       # Standard output
#SBATCH --error=/home/mfafouti/scratch/Mommybrain_marlen/Slide_tags/Single_cell_analysis/Integration/Integration_out/%x_%j.err        # Standard error
#SBATCH --mail-user=mariaeleni.fafouti@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

source /home/mfafouti/miniforge3/bin/activate
conda init 
conda activate seurat_env
cd /home/mfafouti/scratch/Mommybrain_marlen/Mommybrain-PPD/Slide_tags
Rscript integration_seurat.R
conda deactivate 
