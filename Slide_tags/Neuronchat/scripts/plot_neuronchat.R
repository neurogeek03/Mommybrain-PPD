# APPTAINER 
# apptainer shell \
#     --bind ./data:/mnt/data \
#     --bind ./out:/mnt/out \
#     --bind ./scripts:/mnt/scripts \
#     /scratch/mfafouti/docker/neuronchat_full.sif

# ================ ARGS ================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript my_script.R <name> <age>")}

sample <- args[1]


library(NeuronChat)
library(CellChat)
library(glue)
library(ggalluvial)
library(optparse)

print('Libraries loaded!')

# ================ ARGS ================
option_list <- list(
    make_option(c("-d", "--data_dir"), type="character", default="/mnt/data/BC13_TEST",
                help="Directory containing the input data", metavar="character"),
    make_option(c("-o", "--output_dir"), type="character", default="/mnt/out",
                help="Directory for output files", metavar="character")
)

# ================ PARAMS ================
# general
class_column  <-"class_name"
subclass_column <- "subclass_name"

#input
obj_dir <- "/mnt/out"
data_dir <- "/mnt/data/BC13_TEST"
meta_path <- file.path(data_dir, "pos_BC13_metadata.csv")
#output
out_dir <- "/mnt/out/figures"
x_name <- "ran_nc_object.rds"
csv_out_name <- "net_aggregated_x.csv"


csv_out_path <- file.path(obj_dir, csv_out_name)
x_rds_filepath <- file.path(obj_dir, x_name)

# Loading ran NC object 
x <- readRDS(x_rds_filepath)
net_aggregated_x <- net_aggregation(x@net,method = 'weight')

# ================ PREPROCESSING ================
# meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
# print(glue("Metadata loaded. Dimensions: {nrow(meta_data)} rows, {ncol(meta_data)} columns."))

# Creating a mapping of subclass names to borader class names
meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
df_group <- meta_data[!duplicated(meta_data[[subclass_column]]), c(class_column, subclass_column)]
group <- structure(df_group[[class_column]], names = df_group[[subclass_column]])

# ================ HEATMAP ================
heatmap_path = file.path(out_dir, "heatmap_aggregated.png")
png(heatmap_path, width = 30, height = 25, units = "in", res = 300)
heatmap_aggregated(x, method='weight',group=group)
dev.off()

# ================ CIRCLES ================
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

# Set the path for the output PDF
rank_plot_path <- file.path(out_dir, "ranked_signaling_plots.pdf")

# Open a PDF device to save the plots
pdf(rank_plot_path, width = 8, height = 10)

# Plot 1: Rank OUTGOING signaling strength for each pathway
# 'stacked = FALSE' shows the total strength from all cell groups
rankNet_Neuron(x,
    mode = "single",
    measure = "weight",
    stacked = FALSE,
    do.flip = TRUE,
    title = "Outgoing Signaling Strength (All Pathways)"
)

# Plot 2: Rank INCOMING signaling strength for each pathway
# 'stacked = TRUE' shows the contribution of each cell group to the pathway's strength
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

