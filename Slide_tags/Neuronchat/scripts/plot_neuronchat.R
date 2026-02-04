# APPTAINER 
# apptainer shell \
#     --bind ./data:/mnt/data \
#     --bind ./out:/mnt/out \
#     --bind ./scripts:/mnt/scripts \
#     /scratch/mfafouti/docker/neuronchat_full.sif
library(NeuronChat)
library(CellChat)
library(glue)
library(ggalluvial)
# library(optparse)

print('Libraries loaded!')

# ================ ARGS ================
# This script is designed to be called from the command line.
# It expects two arguments:
# --data_dir: The directory containing the input files for a specific sample.
#             This directory should have the ran_nc_object.rds and the metadata file.
# --sample: The sample ID (e.g., "BC13").

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
    stop("Usage: Rscript plot_neuronchat.R --data_dir <path_to_sample_data> --sample <sample_id>", call. = FALSE)
}

# ================ PATHS AND PARAMS ================
# General params
class_column  <- "class_name"
subclass_column <- "subclass_name"

# Input paths
# The ran neuronchat object is expected to be in the data_dir
# Note: The 'test_neuronchat.R' script saves '..._neuronchat_object.rds'. 
# This plotting script assumes you have run the next step and saved the result as 'ran_nc_object.rds'.
# If your file is named differently, you may need to adjust 'x_name'.
x_name <- glue("{sample}_neuronchat_object.rds")
x_rds_filepath <- file.path(data_dir, x_name)

# Metadata is also expected in the data_dir
meta_name <- glue("{sample}_metadata.csv")
meta_path <- file.path(data_dir, meta_name)

# Output paths
# Plots will be saved in a 'figures' subdirectory inside the data_dir
out_dir <- file.path(data_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

print(glue("Input object: {x_rds_filepath}"))
print(glue("Metadata file: {meta_path}"))
print(glue("Output directory for plots: {out_dir}"))


# ================ DATA LOADING ================
# Loading ran NC object 
if (!file.exists(x_rds_filepath)) {
    stop(glue("ERROR: Input file not found at {x_rds_filepath}"), call. = FALSE)
}
x <- readRDS(x_rds_filepath)

net_aggregated_x <- net_aggregation(x@net,method = 'weight')

# ================ PREPROCESSING ================
# meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
# print(glue("Metadata loaded. Dimensions: {nrow(meta_data)} rows, {ncol(meta_data)} columns."))

# Loading metadata to create cell type groupings
if (!file.exists(meta_path)) {
    stop(glue("ERROR: Metadata file not found at {meta_path}"), call. = FALSE)
}
meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
df_group <- meta_data[!duplicated(meta_data[[subclass_column]]), c(class_column, subclass_column)]
group <- structure(df_group[[class_column]], names = df_group[[subclass_column]])

# Aggregated network for circle plot
net_aggregated_x <- net_aggregation(x@net, method = 'weight')

# ================ PLOTTING: HEATMAP ================
print("Generating heatmap...")
heatmap_path = file.path(out_dir, "heatmap_aggregated.png")
png(heatmap_path, width = 30, height = 25, units = "in", res = 300)
heatmap_aggregated(x, method='weight', group=group)
dev.off()
print(glue("  - Saved to {heatmap_path}"))


# ================ PLOTTING: CIRCLE PLOT ================
print("Generating circle plot...")
top = 0.05
circles_path = file.path(out_dir, glue("arrow_circle_nc_{top}.png"))
png(circles_path, width = 30, height = 25, units = "in", res = 300)
# par(mfrow=c(1,2))
# Visualization, circle plot, for the aggregated network
netVisual_circle_neuron(net_aggregated_x,group=group,vertex.label.cex = 1, top = top)
# Visualization, chordDiagram, for the aggregated network; also using cellchat function netVisual_chord_cell_internal(net_aggregated_x, group = group,lab.cex=1)
# CellChat::netVisual_chord_cell_internal(net_aggregated_x, group = group, link.lwd = 2)
# netVisual_chord_neuron(x,method = 'weight',group=group,lab.cex = 1.5)
dev.off()

# ================ RANK INTERACTIONS ================
# This generates plots ranking the signaling pathways by their overall strength.

# ================ PLOTTING: RANKED INTERACTIONS ================
print("Generating ranked signaling plots...")
rank_plot_path <- file.path(out_dir, "ranked_signaling_plots.pdf")

# Open a PDF device to save the plots
pdf(rank_plot_path, width = 8, height = 10)

# Plot 1: Rank OUTGOING signaling strength
rankNet_Neuron(x,
    mode = "single",
    measure = "weight",
    stacked = FALSE,
    do.flip = TRUE,
    title = "Outgoing Signaling Strength (All Pathways)"
)

# Plot 2: Rank INCOMING signaling strength
rankNet_Neuron(x,
    mode = "single",
    measure = "weight",
    stacked = TRUE,
    do.flip = TRUE,
    title = "Incoming Signaling Strength (All Pathways)"
)

# Close the PDF device
dev.off()


# print("--- Inspecting the 'net_aggregated_x' object ---")
# cat("Class of object:", class(net_aggregated_x), "\n")
# cat("\nStructure of object (str):\n")
# str(net_aggregated_x)
# cat("\n")
# print("--- done ---")


# output_png_file <- "/mnt/out/heatmap_aggregated.png"
# png(filename = output_png_file, width = 800, height = 700, res = 100)
# heatmap_aggregated(x, method='weight',group=group, cut_off = 0.05)
# dev.off()

# output_png_file <- "/mnt/out/net_aggregated_x.png"

# # cell_subtypes_to_keep <- c(
# #     "334 Microglia NN",
# #     "319 Astro-TE NN",
# #     "016 CA1-ProS Glut",
# #     "017 CA3 Glut", 
# #     "027 L6b EPd Glut")

# filtered_group <- group[names(group) %in% cell_subtypes_to_keep]
# print(filtered_group)
# netVisual_circle_neuron(net_aggregated_x,group=group,vertex.label.cex = 0.4, top = 0.005, remove.isolate = TRUE)
# # Visualization, chordDiagram, for the aggregated network; also using cellchat function netVisual_chord_cell_internal(net_aggregated_x, group = group,lab.cex=1)

