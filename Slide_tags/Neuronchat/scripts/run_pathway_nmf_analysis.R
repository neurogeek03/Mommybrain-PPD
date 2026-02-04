# ================ PACKAGES ================ 
library(NeuronChat)
library(NMF)
library(ComplexHeatmap)
library(ggplot2)
library(glue)
library(methods) 
library(circlize)
library(ggalluvial) # Required for river plots

print('Libraries loaded for pathway NMF analysis!')

# ================ ARGS ================ 
args <- commandArgs(trailingOnly = TRUE)

input_neuronchat_object_path <- NULL
output_dir <- NULL
k_nmf <- NULL
k_range_str <- NULL

i <- 1
while (i <= length(args)) {
    if (args[i] == "--input_neuronchat_object_path" && i + 1 <= length(args)) {
        input_neuronchat_object_path <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--output_dir" && i + 1 <= length(args)) {
        output_dir <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--k" && i + 1 <= length(args)) {
        k_nmf <- as.numeric(args[i + 1])
        i <- i + 2
    } else if (args[i] == "--k_range" && i + 1 <= length(args)) {
        k_range_str <- args[i + 1]
        i <- i + 2
    } else {
        i <- i + 1
    }
}

if (is.null(input_neuronchat_object_path) || is.null(output_dir)) {
  stop("Usage: Rscript run_pathway_nmf_analysis.R --input_neuronchat_object_path <path> --output_dir <path> [--k <number> | --k_range <start:end>]", call. = FALSE)
}
if (is.null(k_nmf) && is.null(k_range_str)) {
  stop("You must provide either --k (to run NMF) or --k_range (to estimate the optimal k).", call. = FALSE)
}
if (!is.null(k_nmf) && !is.null(k_range_str)) {
  stop("You cannot provide both --k and --k_range. Please choose one.", call. = FALSE)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ================ LOAD AND PREPARE DATA ================ 
print(glue("Loading NeuronChat object from {input_neuronchat_object_path}"))
neuronchat_object <- readRDS(input_neuronchat_object_path)

print("Preparing pathway-interaction matrix...")
net_list <- methods::slot(neuronchat_object, "net")

all_cell_groups <- sort(unique(unlist(lapply(net_list, function(mat) c(rownames(mat), colnames(mat))))))

pathway_interaction_matrix_list <- lapply(names(net_list), function(pathway_name) {
  mat <- net_list[[pathway_name]]
  full_mat <- matrix(0, nrow = length(all_cell_groups), ncol = length(all_cell_groups),
                     dimnames = list(all_cell_groups, all_cell_groups))
  
  if (nrow(mat) > 0 && ncol(mat) > 0) {
      full_mat[rownames(mat), colnames(mat)] <- mat
  }

  as.vector(full_mat)
})

pathway_interaction_matrix_raw <- do.call(rbind, pathway_interaction_matrix_list)
rownames(pathway_interaction_matrix_raw) <- names(net_list)

zero_rows <- which(rowSums(pathway_interaction_matrix_raw) == 0)
if (length(zero_rows) > 0) {
  pathway_interaction_matrix <- pathway_interaction_matrix_raw[-zero_rows, , drop = FALSE]
} else {
  pathway_interaction_matrix <- pathway_interaction_matrix_raw
}

zero_cols <- which(colSums(pathway_interaction_matrix) == 0)
if (length(zero_cols) > 0) {
  pathway_interaction_matrix <- pathway_interaction_matrix[, -zero_cols, drop = FALSE]
}


# =======================================================
#               --- CHOOSE SCRIPT MODE --- 
# =======================================================

if (!is.null(k_range_str)) {
    # --- MODE 1: SELECT K --- 
    print("--- Running in 'select k' mode ---")
    
    k_range_parts <- as.numeric(strsplit(k_range_str, ":")[[1]])
    k_range <- k_range_parts[1]:k_range_parts[2]

    print(glue("Estimating rank for k in range {min(k_range)} to {max(k_range)}..."))

    res <- NMF::nmfEstimateRank(pathway_interaction_matrix, range = k_range, method = 'lee', nrun = 30, seed = 123)

    df1 <- data.frame(k = res$measures$rank, score = res$measures$cophenetic, Measure = "Cophenetic")
    df2 <- data.frame(k = res$measures$rank, score = res$measures$silhouette.consensus, Measure = "Silhouette")
    df <- rbind(df1, df2)

    gg <- ggplot(df, aes(x = k, y = score, group = Measure, color = Measure)) + 
      geom_line(linewidth = 1) + 
      geom_point() + 
      theme_classic() + 
      labs(x = 'Number of patterns (k)', y = 'Measure score', title = "NMF Rank Estimation for Pathways") + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right") + 
      scale_x_continuous(breaks = k_range) + 
      facet_wrap(~ Measure, scales = 'free_y')

    plot_path <- file.path(output_dir, "select_k_plot.pdf")
    ggsave(plot_path, plot = gg, width = 8, height = 5)

    print(glue("Saved plot to {plot_path}"))
    print("---")
    print("Please inspect the plot to select the optimal 'k'. A suitable 'k' is often where the Cophenetic and Silhouette scores begin to drop.")
    print("Once you have chosen a value for 'k', re-run this script using the '--k <your_chosen_k>' argument.")
    print("---")

} else if (!is.null(k_nmf)) {
    # --- MODE 2: RUN NMF ANALYSIS --- 
    print(glue("--- Running in 'NMF analysis' mode with k = {k_nmf} ---"))
    
    if (nrow(pathway_interaction_matrix) < k_nmf) {
      stop(glue("Number of pathways ({nrow(pathway_interaction_matrix)}) is less than k ({k_nmf}). Cannot perform NMF. Adjust k or provide more data."), call. = FALSE)
    }
    if (ncol(pathway_interaction_matrix) < k_nmf) {
      stop(glue("Number of interaction features ({ncol(pathway_interaction_matrix)}) is less than k ({k_nmf}). Cannot perform NMF. Adjust k or provide more data."), call. = FALSE)
    }
    if (k_nmf <= 1) {
      stop("k must be greater than 1 for NMF.", call. = FALSE)
    }

    # --- Custom Pathway NMF Analysis --- 
    print("--- Running Custom Pathway NMF Analysis ---")
    set.seed(123) # for reproducibility
    nmf_res <- NMF::nmf(pathway_interaction_matrix, rank = k_nmf, method = "lee", seed = "nndsvd")

    W <- NMF::basis(nmf_res) # Pathways x Patterns
    H <- NMF::coef(nmf_res)  # Patterns x Interactions
    
    colnames(W) <- paste0("Pattern_", 1:ncol(W))

    # Plotting Logic 1: Pathway Contributions (W matrix)
    print("Generating heatmap for pathway contributions to patterns...")
    col_fun_w = circlize::colorRamp2(c(min(W), median(W), max(W)), c("blue", "white", "red"))
    heatmap_file_path <- file.path(output_dir, glue("pathway_patterns_heatmap_k{k_nmf}.pdf"))
    pdf(heatmap_file_path, width = 8, height = max(4, 0.15 * nrow(W) + 2))
    ht_w <- ComplexHeatmap::Heatmap(W, name = "Contribution",
                                  col = col_fun_w,
                                  cluster_rows = TRUE, cluster_columns = FALSE,
                                  show_row_dend = TRUE, show_column_dend = FALSE,
                                  row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 10),
                                  row_title = "Signaling Pathways", column_title = glue("Pathway NMF Patterns (k={k_nmf})"),
                                  heatmap_legend_param = list(title = "NMF Weight", direction = "horizontal"))
    draw(ht_w, heatmap_legend_side = "bottom")
    dev.off()
    print(glue("Saved Pathway Patterns Heatmap to {heatmap_file_path}"))

    # Plotting Logic 2: Interaction Contributions (H matrix)
    print("Generating heatmaps for cell-cell interaction patterns...")
    H_full <- matrix(0, nrow = k_nmf, ncol = length(all_cell_groups)^2)
    if (length(zero_cols) > 0) {
        H_full[, -zero_cols] <- H 
    } else {
        H_full <- H
    }
    interaction_heatmap_path <- file.path(output_dir, glue("interaction_patterns_heatmap_k{k_nmf}.pdf"))
    pdf(interaction_heatmap_path, width = 8, height = 7)
    for (i in 1:k_nmf) {
        interaction_matrix <- matrix(H_full[i, ], nrow = length(all_cell_groups), ncol = length(all_cell_groups),
                                     dimnames = list(all_cell_groups, all_cell_groups), byrow = FALSE)
        col_fun_h = circlize::colorRamp2(c(min(interaction_matrix), median(interaction_matrix), max(interaction_matrix)), c("ivory", "peachpuff", "red"))
        ht_h <- ComplexHeatmap::Heatmap(interaction_matrix, name = "Contribution", col = col_fun_h,
                                        cluster_rows = TRUE, cluster_columns = TRUE,
                                        row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 8),
                                        row_title = "Sender", column_title = "Receiver", column_title_side = "bottom",
                                        heatmap_legend_param = list(title = "Interaction Weight"))
        draw(ht_h, column_title = glue("Cell-Cell Interaction Profile for Pattern {i}"), column_title_gp = grid::gpar(fontsize = 12))
    }
    dev.off()
    print(glue("Saved Interaction Pattern Heatmaps to {interaction_heatmap_path}"))


    # --- Standard NeuronChat NMF and River Plot Analysis --- 
    print("--- Running Standard NeuronChat NMF for River Plots ---")
    
    # Outgoing patterns
    print("Identifying and plotting 'outgoing' communication patterns...")
    neuronchat_object <- identifyCommunicationPatterns_Neuron(neuronchat_object, pattern = "outgoing", k = k_nmf, heatmap.show = FALSE)
    outgoing_river_plot <- netAnalysis_river_Neuron(neuronchat_object, pattern = "outgoing")
    outgoing_river_plot_path <- file.path(output_dir, glue("outgoing_river_plot_k{k_nmf}.pdf"))
    ggsave(outgoing_river_plot_path, plot = outgoing_river_plot, width = 12, height = 8)
    print(glue("Saved outgoing river plot to {outgoing_river_plot_path}"))

    # Incoming patterns
    print("Identifying and plotting 'incoming' communication patterns...")
    neuronchat_object <- identifyCommunicationPatterns_Neuron(neuronchat_object, pattern = "incoming", k = k_nmf, heatmap.show = FALSE)
    incoming_river_plot <- netAnalysis_river_Neuron(neuronchat_object, pattern = "incoming")
    incoming_river_plot_path <- file.path(output_dir, glue("incoming_river_plot_k{k_nmf}.pdf"))
    ggsave(incoming_river_plot_path, plot = incoming_river_plot, width = 12, height = 8)
    print(glue("Saved incoming river plot to {incoming_river_plot_path}"))

}

print("\n--- Script complete! ---")
