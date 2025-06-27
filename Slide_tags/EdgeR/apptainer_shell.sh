module load apptainer 

apptainer shell \
  --bind /gpfs/fs0/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags/EdgeR:/workspace \
  --bind /gpfs/fs0/scratch/s/shreejoy/mfafouti/miniforge3/envs/rctd_env:/opt/rctd_env \
  --bind /gpfs/fs0/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_seq/Reference_data/new:/project \
  /gpfs/fs0/scratch/s/shreejoy/mfafouti/Mommybrain/edger.sif
