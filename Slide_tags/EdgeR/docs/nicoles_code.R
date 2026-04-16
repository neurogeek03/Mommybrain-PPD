# run on conda activate de_env
setwd("P1_SCZ_DE")
library(ggrepel)
library(cowplot)
library(limma)
library(dplyr)
library(edgeR)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)

meta <- read.csv("Files/Pseudobulk_metadata_supertypes_7cohorts.csv")
bulk <- readRDS("Files/MSSM_pseudobulk_SCZ_supertype.rds")


meta <- meta %>%
  mutate(Donor = gsub("_", "-", Donor)) %>%
  filter(Cohort == "MSSM")  # Keep only MSSM donors

# Get only columns starting with "Sst-" (uppercase S)
sst_cols <- grep("^Sst-", colnames(bulk), value = TRUE)

# Extract unique subtypes (including the number, e.g., "Sst-25")
supertypes <- unique(sub("_.*", "", sst_cols))

print(supertypes)

# Loop over each supertype
for (stype in supertypes) {
  
  message("Processing supertype: ", stype)

  # Identify matching column in metadata (Sst_*)
  count_col <- gsub("-", "_", stype)

  # Keep donors with >=5 cells of this subtype
  donors_keep <- meta %>%
    filter(.data[[count_col]] >= 5) %>%
    pull(Donor)

  # Subset bulk for current supertype
  bulk_subset <- bulk[, grep(paste0("^", stype, "_"), colnames(bulk))]
  
  if (ncol(bulk_subset) == 0) next
  
  pseudobulk <- bulk_subset %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID_Celltype") %>%
    separate(ID_Celltype, into = c("supertype", "Donor"), sep = "_", extra = "merge", fill = "right") %>%
    relocate(supertype, Donor) %>%
    mutate(Donor = gsub("_", "-", Donor)) %>%
    filter(Donor %in% donors_keep)
  
  # Create matrix for DGEList
  mat <- pseudobulk %>%
    tibble::column_to_rownames("Donor") %>%
    dplyr::select(-supertype) %>%
    as.matrix() %>%
    t()
  
  # Reorder metadata to match columns
  meta_subset <- meta[match(colnames(mat), meta$Donor), ]
  
     # Skip if Diagnosis has <2 levels
  if (length(unique(meta_subset$Diagnosis)) < 2) {
    message("Skipping ", stype, " because Diagnosis has <2 levels")
    next
  }
  
  # DGEList
  dge <- DGEList(counts = mat, genes = rownames(mat))
  
  # Keep genes with ≥1 count in ≥80% of samples
  min_samples <- ncol(mat) * 0.8
  dge <- dge[rowSums(dge$counts >= 1) >= min_samples, ]
  
  # Normalization
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Design matrix
  design <- model.matrix(~ scale(Age) + Sex + Diagnosis + scale(PMI), data = meta_subset)
  
  # voom + lmFit + eBayes
  vm <- voom(dge, design, plot = FALSE)
  fit <- lmFit(vm, design)
  fit <- eBayes(fit)
  
  # Extract DE results
  DE <- topTable(
    fit,
    coef = "DiagnosisSchizophrenia",
    n = Inf,
    adjust.method = "BH"
  )
  
  
  # Save results per supertype
  saveRDS(DE, paste0("Files/DE_results_MSSM_", stype, ".rds"))
  
}