# =====================================================================
# Title:        Doublet detection 
# Description:  The integrated seurat object is processed to remove 
# doublets, using a custom function with the DoubletFinder package
# Author:       Maria Eleni Fafouti
# Date:         2025-04-28
# =====================================================================

# ---- 1. SETUP ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(tibble)

set.seed(27)

# Define your multiplet rate here. For my research, I consulted the table from: https://cdn.10xgenomics.com/image/upload/v1710230393/support-documents/CG000731_ChromiumGEM-X_SingleCell3_ReagentKits_v4_UserGuide_RevA.pdf
percent_of_multiplets = 0.05

# Relevant paths
project_path <- "/scratch/s/shreejoy/mfafouti/Mommybrain/Slide_tags"
in_path <- file.path(project_path, 'Integration', 'Integration_CB_out') 
doublet_folder <- file.path(project_path, 'Doublet_detection', 'dblfinder_1st_try') 

# Loading the function
source(file.path(doublet_folder, "run_doubletfinder_custom.R"))

# Read in your Seurat object
seu <- readRDS(file.path(in_path, "new_integrated.rds"))
message("Seurat object loaded!")


# ---- 2. PRE-PROCESSING THE SEURAT OBJECT ----
# Previewing metadata
head(seu@meta.data)
table(seu$orig.ident)

# Joining layers 
seu[["RNA"]] <- JoinLayers(seu[["RNA"]]) 
message("Layers Joined!")

# Switching to the RNA assay (raw data) - not the integrated one - because DoubletFinder requires unintegrated data
DefaultAssay(seu) <- "RNA"
message("Assay switched!")


# ---- 3. RUNNING DOUBLET FINDER ----
# DoubletFinder should be run on a per sample basis, so we split the object into the individual samples 
samp_split <- SplitObject(seu, split.by = "orig.ident") 
head(samp_split) #just to verify

# Get Doublet/Singlet IDs by DoubletFinder()
samp_split <- lapply(samp_split, run_doubletfinder_custom) 

# Setting the multiplet rate manually:
samp_split <-lapply(samp_split, run_doubletfinder_custom, multiplet_rate = percent_of_multiplets) # here we run doubletfinder 
message("DoubletFinder completed!")
head(samp_split)

# ---- 4. ADDING DOUBLET INFO TO METADATA ----
# merge to a single dataframe
sglt_dblt_metadata <- data.frame(bind_rows(samp_split)) 

# assign cell IDs to row names to ensure match
rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names 
sglt_dblt_metadata$row_names <- NULL
head(sglt_dblt_metadata)

# Adding metadata to seurat object
seu <- AddMetaData(seu, sglt_dblt_metadata, col.name = 'doublet_finder')
head(seu@meta.data)


# ---- 5. ADDING DOUBLET INFO TO METADATA ----

# Violin Plot: Check how doublets singlets differ in QC measures per sample.
# Create a long-format data frame for ggplot
plot_data <- meta %>%
  select(orig.ident, scDblFinder.class, nFeature_RNA, nCount_RNA) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA),
               names_to = "feature", values_to = "value")

vln_plot <- ggplot(plot_data, aes(x = orig.ident, y = value, fill = scDblFinder.class)) +
  geom_violin(scale = "width", trim = FALSE) +
  facet_wrap(~ feature, ncol = 2, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "Value", fill = "Doublet Class")

ggsave(
  filename = file.path(doublet_folder, "out", "scdbl_qc_vlnplot.png"),
  plot = vln_plot,
  width = 20,
  height = 8,
  dpi = 300
)

# Table: Get doublets per sample
meta <- meta %>%
  rename(doublet_finder = scDblFinder.class)

# Compute summary
doublets_summary <- meta %>% 
  group_by(orig.ident, doublet_finder) %>% 
  summarise(total_count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(doublet_finder, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count / countT, 2), "%")) %>%
  select(-countT)
write.table(doublets_summary, file = file.path(doublet_folder, '_doubletfinder_doublets_summary.txt'), quote = FALSE, row.names = FALSE, sep = '\t')

# ---- 6. SAVING DATA ----
# Saving the seurat object
saveRDS(seu, file = file.path(doublet_folder, "sobj_doublets_5%.rds"))
