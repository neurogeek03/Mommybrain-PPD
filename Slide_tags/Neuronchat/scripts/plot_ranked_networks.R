# This script ranks the signaling networks and plots them.

library(NeuronChat)
library(CellChat)
library(glue)
library(RColorBrewer) # Needed if you want to use RColorBrewer palettes

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
# Creating a mapping of subclass names to borader class names
meta_data <- read.csv(meta_path, row.names = 1, check.names = FALSE)
df_group <- meta_data[!duplicated(meta_data[[subclass_column]]), c(class_column, subclass_column)]
group <- structure(df_group[[class_column]], names = df_group[[subclass_column]])

# ================ RANKING NETWORKS ================
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
