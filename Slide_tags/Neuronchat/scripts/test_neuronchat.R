APPTAINER 
apptainer shell \
    --bind ./data:/mnt/data \
    --bind ./out:/mnt/out \
    --bind ./scripts:/mnt/scripts \
    /scratch/mfafouti/docker/neuronchat_full.sif

library(NeuronChat)
library(CellChat)
library(glue)
library(parallel)

print('Libraries loaded!')

# ================ INPUT ================
data_dir <- "/mnt/data/BC13_TEST"
output_dir <- "/mnt/out"
matrix_path <- file.path(data_dir, "BC13_expression_matrix_cell_subclass.csv")
meta_path <- file.path(data_dir, "BC13_metadata.csv")

# ================ PARAMS ================
class_column  <-"class_name"
subclass_column <- "subclass_name"
x_name <- "neuronchat_object.rds"
x_rds_filepath <- file.path(output_dir, x_name)

# # Checking compute resources available
# num_cores <- parallel::detectCores()
# cat("Number of CPU cores detected:", num_cores, "\n")

# # ================ Data loading ================
# expr_matrix <- read.csv(matrix_path, row.names = 1, check.names = FALSE)

# # Ensuring that genes are rows and cells are columns 
# print("Transposing the expression matrix to be [genes x cells]...")
# expr_matrix_transposed <- t(expr_matrix)
# print(glue("Expression matrix loaded. Dimensions: {nrow(expr_matrix)} rows, {ncol(expr_matrix)} columns."))
# meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
# print(glue("Metadata loaded. Dimensions: {nrow(meta_data)} rows, {ncol(meta_data)} columns."))

# # Creating a mapping of subclass names to borader class names 
# df_group <- meta_data[!duplicated(meta_data$subclass_label), c(class_column, subclass_column)]
# group <- structure(df_group$class_label,names=df_group$subclass_label)

# # ================ Create object & run ================
# x <- createNeuronChat(t(as.matrix(expr_matrix[,1:(dim(expr_matrix)[2]-1)])),DB='mouse',group.by = expr_matrix$cell_subclass);
# saveRDS(x, file = x_rds_filepath)

neuronchat_obj_path <- file.path({output_dir}, {x_name})
nc_obj <- readRDS(neuronchat_obj_path)

x <- run_NeuronChat(x,M=100)
# inspect object
slotNames(x)
head(x@meta)
levels(x@idents)
str(x@net)

net_aggregated_x <- net_aggregation(x@net,method = 'weight')

print("--- Inspecting the 'net_aggregated_x' object ---")
cat("Class of object:", class(net_aggregated_x), "\n")
cat("\nStructure of object (str):\n")
str(net_aggregated_x)
cat("\n")
print("--- done ---")