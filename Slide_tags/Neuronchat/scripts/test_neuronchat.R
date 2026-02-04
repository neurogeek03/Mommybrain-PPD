# APPTAINER 
# apptainer shell \
#     --bind ./data:/mnt/data \
#     --bind ./out:/mnt/out \
#     --bind ./scripts:/mnt/scripts \
#     /scratch/mfafouti/docker/neuronchat_full.sif

# ================ PACKAGES ================
library(NeuronChat)
library(CellChat)
library(glue)
library(parallel)

print('Libraries loaded!')

# ================ ARGS ================
args <- commandArgs(trailingOnly = TRUE)

data_dir <- NULL
sample <- NULL

i <- 1
while (i <= length(args)) {
    if (args[i] == "--data_dir" && i + 1 <= length(args)) {
        data_dir <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--sample" && i + 1 <= length(args)) {
        sample <- args[i + 1]
        i <- i + 2
    } else {
        i <- i + 1
    }
}

if (is.null(data_dir) || is.null(sample)) {
    stop("Both --data_dir and --sample arguments must be supplied. Usage: Rscript test_neuronchat.R --data_dir <path> --sample <id>", call. = FALSE)
}

# ================ INPUT ================
matrix_path <- file.path(data_dir, glue("{sample}_expression_matrix_cell_subclass.csv"))
meta_path <- file.path(data_dir,  glue("{sample}_metadata.csv"))

# ================ PARAMS ================
class_column  <-"class_name"
subclass_column <- "subclass_name"
x_name <- glue("{sample}_neuronchat_object.rds")
x_rds_filepath <- file.path(data_dir, x_name)

# Checking compute resources available
num_cores <- parallel::detectCores()
cat("Number of CPU cores detected:", num_cores, "\n")

# ================ Data loading ================
expr_matrix <- read.csv(matrix_path, row.names = 1, check.names = FALSE)
print(glue("Expression matrix loaded. Dimensions: {nrow(expr_matrix)} rows, {ncol(expr_matrix)} columns."))
meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
print(glue("Metadata loaded. Dimensions: {nrow(meta_data)} rows, {ncol(meta_data)} columns."))

# Creating a mapping of subclass names to borader class names 
df_group <- meta_data[!duplicated(meta_data$subclass_label), c(class_column, subclass_column)]
group <- structure(df_group$class_name, names=df_group$subclass_name)
# group <- structure(df_group$class_label,names=df_group$subclass_label)

# # ================ Create object & run ================
# x <- createNeuronChat(t(as.matrix(expr_matrix[,1:(dim(expr_matrix)[2]-1)])),DB='mouse',group.by = expr_matrix$cell_subclass);
# saveRDS(x, file = x_rds_filepath)

# ================ Data loading ================
expr_df <- read.csv(matrix_path, row.names = 1, check.names = FALSE)
print(glue("Loaded data with dimensions: {nrow(expr_df)} rows, {ncol(expr_df)} columns."))

# Ensuring that genes are rows and cells are columns 
# print("Transposing the expression matrix to be [genes x cells]...")
# expr_matrix <- t(expr_matrix)
# print(glue("Expression matrix loaded. Dimensions: {nrow(expr_matrix)} rows, {ncol(expr_matrix)} columns."))
cell_subclass_labels <- expr_df$cell_subclass
expr_matrix <- as.matrix(expr_df[, -ncol(expr_df)]) # numeric matrix [cells x genes]
print("Separated cell subclass labels from the expression matrix.")
# Transpose the expression matrix to be [genes x cells] for NeuronChat
expr_matrix_t <- t(expr_matrix)
print(glue("Transposed expression matrix to dimensions: {nrow(expr_matrix_t)} rows, {ncol(expr_matrix_t)} columns."))

meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
print(glue("Metadata loaded. Dimensions: {nrow(meta_data)} rows, {ncol(meta_data)} columns."))

# Creating a mapping of subclass names to borader class names 
df_group <- meta_data[!duplicated(meta_data$subclass_label), c(class_column, subclass_column)]
group <- structure(df_group$class_label,names=df_group$subclass_label)

# ================ Create object & run ================
x <- createNeuronChat(expr_matrix_t, meta=meta_data, DB='mouse', group.by = cell_subclass_labels);
print("Running main NeuronChat analysis (this may take a while)...")

x <- run_NeuronChat(x, M=100) 
print(glue("Analysis complete. Saving processed object to {x_rds_filepath}"))

saveRDS(x, file = x_rds_filepath)

# # ================ BEGIN FROM SAVED OBJECT ================
# neuronchat_obj_path <- file.path({output_dir}, 'auto','test5', {x_name})
# x <- readRDS(neuronchat_obj_path)
# # ================

# x <- run_NeuronChat(x,M=100)

# ran_nc_filepath = file.path(output_dir, 'ran_nc_object.rds')
# saveRDS(x, file = ran_nc_filepath)

# # inspect object
# slotNames(x)
# head(x@meta)
# levels(x@idents)
# str(x@net)

