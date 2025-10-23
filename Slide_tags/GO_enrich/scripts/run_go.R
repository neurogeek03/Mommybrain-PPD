library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(cowplot)
library(ggplot2)
library(dplyr) 
library(tidyr)
library(stringr)

# de <- readRDS("Files/DE_results_Affected_vs_Unaffected.rds")
# sig_de <- subset(de, adj.P.Val < 0.05)

# # Go enrichment
# ego <- enrichGO(gene          = sig_de$genes,
#                    OrgDb         = org.Hs.eg.db,
#                    keyType       = "SYMBOL",
#                    ont           = "BP",          # Biological Process
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.05,
#                    qvalueCutoff  = 0.05)



# p <- barplot(ego, showCategory = 15, title = "GO: Affected v Unaffected")
# ggsave("Figures/GO_Affected_vs_Unaffected2.png", p, width = 14, height = 10, dpi = 300)

