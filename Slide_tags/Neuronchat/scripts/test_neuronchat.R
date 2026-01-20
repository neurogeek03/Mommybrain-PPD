# APPTAINER 
# apptainer shell \
#     --bind ./data:/mnt/data \
#     --bind ./out:/mnt/out \
#     --bind ./scripts:/mnt/scripts \
#     /scratch/mfafouti/docker/neuronchat_full.sif

library(NeuronChat)
library(CellChat)
library(glue)
library(parallel)

print('Libraries loaded!')

# ================ INPUT ================
data_dir <- "/mnt/data/BC13_TEST"
output_dir <- "/mnt/out"
matrix_path <- file.path(data_dir, "pos_BC13_expression_matrix_cell_subclass.csv")
meta_path <- file.path(data_dir, "pos_BC13_metadata.csv")

# ================ PARAMS ================
class_column  <-"class_name"
subclass_column <- "subclass_name"
x_name <- "neuronchat_object.rds"
x_rds_filepath <- file.path(output_dir, x_name)

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

# ================ Create object & run ================
x <- createNeuronChat(t(as.matrix(expr_matrix[,1:(dim(expr_matrix)[2]-1)])),DB='mouse',group.by = expr_matrix$cell_subclass);
saveRDS(x, file = x_rds_filepath)

# # ================ Data loading ================
# expr_df <- read.csv(matrix_path, row.names = 1, check.names = FALSE)
# print(glue("Loaded data with dimensions: {nrow(expr_df)} rows, {ncol(expr_df)} columns."))

# # Ensuring that genes are rows and cells are columns 
# # print("Transposing the expression matrix to be [genes x cells]...")
# # expr_matrix <- t(expr_matrix)
# # print(glue("Expression matrix loaded. Dimensions: {nrow(expr_matrix)} rows, {ncol(expr_matrix)} columns."))
# cell_subclass_labels <- expr_df$cell_subclass
# expr_matrix <- as.matrix(expr_df[, -ncol(expr_df)]) # numeric matrix [cells x genes]
# print("Separated cell subclass labels from the expression matrix.")
# # Transpose the expression matrix to be [genes x cells] for NeuronChat
# expr_matrix_t <- t(expr_matrix)
# print(glue("Transposed expression matrix to dimensions: {nrow(expr_matrix_t)} rows, {ncol(expr_matrix_t)} columns."))

# meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
# print(glue("Metadata loaded. Dimensions: {nrow(meta_data)} rows, {ncol(meta_data)} columns."))

# # Creating a mapping of subclass names to borader class names 
# df_group <- meta_data[!duplicated(meta_data$subclass_label), c(class_column, subclass_column)]
# group <- structure(df_group$class_label,names=df_group$subclass_label)

# # ================ Create object & run ================
# x <- createNeuronChat(expr_matrix_t, DB='mouse', group.by = cell_subclass_labels);
# saveRDS(x, file = x_rds_filepath)

# ================ BEGIN FROM SAVED OBJECT ================
neuronchat_obj_path <- file.path({output_dir}, 'auto','test5', {x_name})
x <- readRDS(neuronchat_obj_path)
# ================

x <- run_NeuronChat(x,M=100)

ran_nc_filepath = file.path(output_dir, 'ran_nc_object.rds')
saveRDS(x, file = ran_nc_filepath)

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


output_png_file <- "/mnt/out/heatmap_aggregated.png"
png(filename = output_png_file, width = 800, height = 700, res = 100)
heatmap_aggregated(x, method='weight',group=group, cut_off = 0.05)
dev.off()

output_png_file <- "/mnt/out/net_aggregated_x.png"

cell_subtypes_to_keep <- c(
    "334 Microglia NN",
    "319 Astro-TE NN",
    "016 CA1-ProS Glut",
    "017 CA3 Glut", 
    "027 L6b EPd Glut")

filtered_group <- group[names(group) %in% cell_subtypes_to_keep]
print(filtered_group)
netVisual_circle_neuron(net_aggregated_x,group=group,vertex.label.cex = 0.4, top = 0.005, remove.isolate = TRUE)
# Visualization, chordDiagram, for the aggregated network; also using cellchat function netVisual_chord_cell_internal(net_aggregated_x, group = group,lab.cex=1)

