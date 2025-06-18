# Mommybrain: Computational genomics on the postpartum brain 
Spatial Transcriptomics (Slide-seq V2, Slide-Tags) project on Postpartum Depression (PPD). Brains of female rats were sequenced at different stages (nulliparous, pregnant, postpartum), while some animals were treated with corticosterone to mimic the psychiatric symptoms of PPD.  

## Explanation of Methods 
![image](https://github.com/user-attachments/assets/be8514d0-ad38-4e93-8ab4-0b7a19070295)

**Slide-seq**: a spatial bead can capture 1 or more cells, collecting information from both the cytoplasm and the nucleus

**Slide-tags**: a spatial bead can capture 1 cell only, collecting information from its nucleus

## Containers used: 
### Seurat v5 container
An easy way to access many Seurat-related packages. Additional packages can be added by mounting a conda env to the container. 
To obtain the Docker image, I run the following on my HPC (which uses Apptainer):
```bash
module load apptainer
apptainer pull seurat_v5.sif docker://satijalab/seurat:5.0.0
```
