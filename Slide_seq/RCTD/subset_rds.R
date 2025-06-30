library(Seurat)

input_dir <- "/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/renamed_rds_files"
output_dir <- "/home/mfafouti/scratch/Mommybrain_marlen/Slide_seq/rctd_run/subsets_100"
num_subsets <- 100  # Replace with however many equal parts you want

# ==== LOAD OBJECTS AND PROCESS ====
files <- list.files(input_dir, pattern = "*.rds", full.names = TRUE)

for (input_file in files) {
  sample <- tools::file_path_sans_ext(basename(input_file))  # Extract sample name, e.g., B08_renamed_seurat → B08_renamed_seurat
  sample_name <- sub("_.*", "", sample)  # Extract 'B08' from 'B08_renamed_seurat'

  seurat_obj <- readRDS(input_file)
  all_cells <- colnames(seurat_obj)

  # ==== SHUFFLE AND SPLIT ====
  set.seed(123)  # For reproducibility
  shuffled_cells <- sample(all_cells)
  split_cells <- split(shuffled_cells, cut(seq_along(shuffled_cells), num_subsets, labels = FALSE))

  # ==== CREATE SAMPLE-SPECIFIC OUTPUT DIR ====
  sample_output_dir <- file.path(output_dir, sample_name)
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)

  # ==== SUBSET & SAVE ====
  for (i in seq_along(split_cells)) {
    sub_obj <- subset(seurat_obj, cells = split_cells[[i]])
    saveRDS(sub_obj, file = file.path(sample_output_dir, paste0(sample_name, "_subset_", i, ".rds")))
  }
}