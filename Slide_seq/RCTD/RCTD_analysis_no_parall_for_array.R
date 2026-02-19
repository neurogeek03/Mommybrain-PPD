
# # To run inside the container

# # # bash
# apptainer shell \
#   --bind /home/mfafouti/miniforge3/envs/rctd_env:/opt/rctd_env \
#   --bind /home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run:/workspace \
#   /home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/seurat_v5.sif
# R --no-init-file

# ========= INPUT ARGUMENTS =========
args <- commandArgs(trailingOnly = TRUE)
seeker_file <- args[1]     # e.g. "subset_1.rds"
ref_path     <- args[2]     # e.g. "ref_slide_tags.rds"
output_file <- args[3]
output_obj <- args[4]
output_RCTD <- args[5]

# ========= LOAD FILES =========
cat("Reading Seeker file:", seeker_file, "\n")
seurat_obj <- readRDS(seeker_file)

cat("Reading Reference:", ref_path, "\n")
reference <- readRDS(ref_path)

Sys.setenv(OPENBLAS_NUM_THREADS = 1)  

# ========= LIBRARIES =========
#NOTE: Updated library calling to be compatible with functional docker image
library(Seurat)
packageVersion("Seurat")
library(spacexr, lib.loc = "/opt/miniforge3/envs/rctd_env/lib/R/library")
library(quadprog, lib.loc = "/opt/miniforge3/envs/rctd_env/lib/R/library")
# library(doParallel,lib.loc = "/opt/rctd_env/lib/R/library")
library(stringr)
library(glue)

# # ========== PARALLELIZATION ==========
# library(foreach, lib.loc = "/opt/rctd_env/lib/R/library")
# library(iterators, lib.loc = "/opt/rctd_env/lib/R/library")
# library(doParallel, lib.loc = "/opt/rctd_env/lib/R/library")

# # Use cores only within one machine
# n_cores <- 8
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)

# # Optional: confirm it's registered
# cat("Registered parallel backend with", getDoParWorkers(), "cores.\n")

# ========== PARAMETERS ==========
sample <- "B08"

# ========= PATHS =========
# seeker_file <- glue("/workspace/{sample}_seurat_renamed.rds")
# ref_path <- file.path("/workspace/ref_slide_tags.rds")

# ========= SET UP REFERENCE =========
ref <- readRDS(ref_path)
# ref <- UpdateSeuratObject(ref) # update seurat obj to new Seurat version format 
Idents(ref) <- "MapMyCells_subclass_name"

# ========= ASSIGNING UNIQUE BARCODES  =========
barcodes <- ref@meta.data$barcode  # or your barcode column

# Check duplicates
sum(duplicated(barcodes))

# Make unique if needed
barcodes_unique <- make.unique(barcodes)

# Assign 
counts <- ref[["RNA"]]$counts
colnames(counts) <- barcodes_unique
colnames(ref) <- barcodes_unique
names(ref$MapMyCells_subclass_name) <- barcodes_unique
names(ref$nCount_RNA) <- barcodes_unique

# Confirm
all(colnames(counts) == names(ref$MapMyCells_subclass_name))
all(colnames(counts) == names(ref$nCount_RNA))
all(colnames(counts) == colnames(ref))

# ========= EXTRACT REF INFO FOR RCTD FUNCTION =========
# extract information to pass to the RCTD Reference function
cluster <- as.factor(ref$MapMyCells_subclass_name)

# cleaning cluster names (to remove forward slashes)
cluster_clean <- gsub("[^A-Za-z0-9_]", "_", cluster)  # replace non-alphanum with '_'
cluster <- as.factor(cluster_clean)
names(cluster) <- colnames(ref)

nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)

# creates an RCTD reference object, which summarizes your single-cell 
# reference data in a form that RCTD uses for spatial deconvolution
reference <- Reference(counts, cluster, nUMI, n_max_cells = 5000) 

# ========= RESHAPE SLIDE-SEQ DATA =========
# TODO this can loop over all samples 
slide.seq <- readRDS(seeker_file)

# ========= QUERY SLIDE-SEQ DATA =========
# counts <- slide.seq[["Spatial"]]$counts => in my case i dont have the spatial data as an assay
counts <- LayerData(slide.seq, assay = "RNA", layer = "counts")
# coords <- GetTissueCoordinates(slide.seq)

# Extract spatial info and making it a df
coords <- Embeddings(slide.seq[["spatial"]])
coords <- as.matrix(Embeddings(slide.seq[["spatial"]]))
colnames(coords) <- c("x", "y")
coords <- as.data.frame(coords)
head(coords)
class(coords)

# sense-check:
all(rownames(coords) == colnames(counts))  # should be TRUE

# coords[is.na(colnames(coords))] <- NULL

query <- SpatialRNA(coords, counts, colSums(counts))

# ========= RUNNING RCTD =========
Sys.setenv(RCTD_NUM_THREADS = "1") # no parallelization
RCTD <- create.RCTD(query, reference, max_cores = 1) # cores was 8 in tutorial but i changed it

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)

# ========= SAVE RCTD RESULTS =========
write.csv(RCTD@results$results_df, file = output_file, row.names = TRUE)
saveRDS(slide.seq, file = output_obj)
saveRDS(RCTD, file = output_RCTD)

# # ========= PLOTTING RCTD ANNOTATIONS =========
# library(patchwork)
# p1 <- SpatialDimPlot(slide.seq, group.by = "first_type") + ggtitle("First Type")
# p2 <- SpatialDimPlot(slide.seq, group.by = "second_type") + ggtitle("Second Type")

# # Save individual plots as PNG with unique file names per subset
# ggsave(paste0("workspace/", seeker_file, "_RCTD_first_type.png"), plot = p1, width = 6, height = 6, dpi = 300)
# ggsave(paste0("workspace/", seeker_file, "_RCTD_second_type.png"), plot = p2, width = 6, height = 6, dpi = 300)

# combined_plot <- p1 | p2
# ggsave(paste0("workspace/", seeker_file, "_RCTD_combined_plot.png"), plot = combined_plot, width = 12, height = 6, dpi = 300)
