# ================ PACKAGES ================ 
library(NeuronChat)
library(ggplot2)
library(patchwork)
library(glue)

print('Libraries loaded for comparative analysis!')

# ================ ARGS ================ 
args <- commandArgs(trailingOnly = TRUE)

base_dir <- NULL
output_dir <- NULL

i <- 1
while (i <= length(args)) {
    if (args[i] == "--base_dir" && i + 1 <= length(args)) {
        base_dir <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--output_dir" && i + 1 <= length(args)) {
        output_dir <- args[i + 1]
        i <- i + 2
    } else {
        i <- i + 1
    }
}

if (is.null(base_dir) || is.null(output_dir)) {
  stop("Usage: Rscript run_comparative_analysis.R --base_dir <path> --output_dir <path>", call. = FALSE)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ================ DEFINE SAMPLES AND GROUPS ================ 
# Define which samples belong to which group
SAMPLES_CONTROL <- c("BC3", "BC9", "BC13")
SAMPLES_TREATMENT <- c("BC14", "BC15", "BC28")
SAMPLES_ALL <- c(SAMPLES_CONTROL, SAMPLES_TREATMENT)

# ================ LOAD AND MERGE (if necessary) ================ 
print("--- Step 1: Loading and Merging Samples ---")
merged_rds_path <- file.path(output_dir, "merged_neuronchat_object.rds")

if (file.exists(merged_rds_path)) {
    print(glue("Found existing merged object at {merged_rds_path}. Loading it."))
    neuronchat.merged <- readRDS(merged_rds_path)
} else {
    print("Merged object not found. Loading individual samples to create it.")
    object.list <- list()
    for (sample_id in SAMPLES_ALL) {
        rds_path <- file.path(base_dir, sample_id, glue("{sample_id}_neuronchat_object.rds"))
        if (!file.exists(rds_path)) {
            stop(glue("ERROR: Input file for sample {sample_id} not found at {rds_path}"), call. = FALSE)
        }
        
        # Assign a descriptive name including the condition
        condition <- if (sample_id %in% SAMPLES_CONTROL) "Control" else "Treatment"
        list_name <- glue("{condition}_{sample_id}")
        
        print(glue("Loading {list_name} from {rds_path}"))
        object.list[[list_name]] <- readRDS(rds_path)
    }
    
    # Merge and save the new object
    print("Merging all samples...")
    neuronchat.merged <- mergeNeuronChat(object.list, add.names = names(object.list))
    
    print(glue("Saving merged object to {merged_rds_path}"))
    saveRDS(neuronchat.merged, file = merged_rds_path)
}

# ================ COMPARATIVE ANALYSIS & PLOTTING ================ 
print("\n--- Step 2: Running Comparative Analyses and Generating Plots ---")

# ----- Analysis 1: Compare Total Interaction Strength ----- 
print("Generating plot for total interaction strength comparison...")
# This compares the first Control vs. the first Treatment sample as an example
# You can change the `group` argument to compare any two samples
group_to_compare <- c(glue("Control_{SAMPLES_CONTROL[1]}"), glue("Treatment_{SAMPLES_TREATMENT[1]}"))
gg_compare_strength <- compareInteractions_Neuron(neuronchat.merged, group = group_to_compare, measure = "weight")
gg_compare_count <- compareInteractions_Neuron(neuronchat.merged, group = group_to_compare, measure = "count")

comparison_plot <- gg_compare_strength + gg_compare_count
comparison_plot_path <- file.path(output_dir, "comparison_total_interactions.pdf")
ggsave(comparison_plot_path, plot = comparison_plot, width = 10, height = 6)
print(glue("  - Saved interaction comparison plot to {comparison_plot_path}"))


# ----- Analysis 2: Functional Similarity Embedding ----- 
print("Generating plot for functional similarity embedding...")
# This visualizes all samples. Samples from the same group should cluster together.
embedding_plot_path <- file.path(output_dir, "comparison_functional_embedding.pdf")
pdf(embedding_plot_path, width = 8, height = 8)
netVisual_embedding_Neuron(neuronchat.merged, type = "functional", title = "Functional Similarity of Signaling")
dev.off()
print(glue("  - Saved embedding plot to {embedding_plot_path}"))


# ----- Analysis 3: Rank Pathways Across All Samples ----- 
print("Generating plot for ranked pathway strength...")
# The stacked bar plot shows the contribution of each sample/group to the pathways.
rank_plot_path <- file.path(output_dir, "comparison_ranked_pathways.pdf")
pdf(rank_plot_path, width = 10, height = 14) # Increased height for better readability
rankNet_Neuron(neuronchat.merged, mode = "merged", stacked = TRUE, do.flip = TRUE, title = "Signaling Pathway Strength Across Samples")
dev.off()
print(glue("  - Saved ranked pathway plot to {rank_plot_path}"))

print("\n--- Comparative analysis script complete! ---")