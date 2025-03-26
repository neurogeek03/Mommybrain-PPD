rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)

# Arguments that will be parsed in the command line when running the script

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
        stop("Sample ID, path to input folder, path to output folder, minimal log10(UMI) are needed as input arguments", call.=FALSE)
} else  {
        sample_id = args[1]
	inputdir = args[2]
	datadir = args[3]
	min_log10_nUMI = args[4] # e.g.1.4
}

if (!dir.exists(datadir)){
	message("-------- Create output dir")
	dir.create(datadir)
}else{
	message("-------- Output dir exist")
}

# PARAMETRIZATION
message("-------- Parameterization")
message("min_log10_UMI=", min_log10_nUMI)
min_nUMI_m <- gsub(".", "pt", min_log10_nUMI, fixed=TRUE)
print(min_nUMI_m)

m=40
n=100
p=5
q=10

# m: Grid size for Step 2.
# n: Grid size for Step 3.
# p: Minimum bead barcodes in a m x m um neighborhood to retain beads.
# q: Minimum bead barcodes in a n x n um neighborhood to retain beads.

message("m=", m)
message("n=", n)
message("p=", p)
message("q=", q)

grid_m=10000/m
grid_n=10000/n


# Step 0: Import seurat object
message("-------- Import seurat object")
seuratO <- readRDS(file.path(inputdir, paste0(sample_id, "_seurat.rds")))


# Step 1: Remove bead barcodes with fewer than min_log10_nUMI log10(UMI)
message("-------- Step 1: Remove bead barcodes with fewer than ", min_log10_nUMI, " log10(UMI)")
seuratO <- subset(seuratO, cells=names(seuratO$log10_nCount_RNA)[seuratO$log10_nCount_RNA >= min_log10_nUMI])


# Step 2: Remove bead barcodes in a m x m (micrometers) neighborhood containing fewer than p bead barcodes
message("-------- Step 2: Remove bead barcodes in a ", m, " x ",  m, " um neighborhood containing fewer than ", p," bead barcodes")
x <- cut(seuratO@reductions$SPATIAL@cell.embeddings[,1], breaks = grid_m, labels = seq(1,grid_m))
y <- cut(seuratO@reductions$SPATIAL@cell.embeddings[,2], breaks = grid_m, labels = seq(1,grid_m))
grouping <- data.frame(rownames(seuratO@reductions$SPATIAL@cell.embeddings), x, y)

#print(head(grouping))
groups <- grouping[,c("x","y")]

#library(dplyr)
density<-dplyr::count(groups %>%
  group_by_all)

print(head(density))

png(file.path(datadir, paste0(sample_id,"_step2_density_",m,"_",p,".png")), width=4.5, height=4.5, res = 300, units = "in")
plot(density(density$n), xlab=paste0("# bead barcodes / (",m,"x",m,") um"), main="")
dev.off()

density_high<-as.data.frame(density[density$n>=p,])
grouping_filtered <- merge(grouping, density_high, by=c("x", "y"))


message("-------- Visualize on/off tissue signals after Step 2 filter")

## ON (#bead barcodes >= p on a mxm um neighborhood) 
png(file.path(datadir, paste0(sample_id,"_step2_UMIgte", min_nUMI_m, "_ontissue_", m, "_", p,".png")), width=7, height=4.5, res = 300, units = "in")
DimPlot(seuratO, reduction = "SPATIAL", cells.highlight = grouping_filtered[,3], cols= "white", raster=FALSE)&
  coord_fixed()&NoAxes()&
  theme(legend.position = "none")&
  ggtitle(paste0("Number of beads: " ,length(grouping_filtered[,3])))
dev.off()

## OFF (#bead barcodes < p on a mxm um neighborhood) 
png(file.path(datadir, paste0(sample_id,"_step2_UMIgte", min_nUMI_m, "_offtissue_", m, "_", p, ".png")), width=7, height=4.5, res = 300, units = "in")
DimPlot(seuratO, reduction = "SPATIAL", cells.highlight = setdiff(grouping[,1], grouping_filtered[,3]), cols= "white", raster=FALSE)&
  coord_fixed()&NoAxes()&
  theme(legend.position = "none")&
  ggtitle(paste0("Number of beads: " ,length(setdiff(grouping[,1], grouping_filtered[,3]))))
dev.off()


# Step 3: Remove bead barcodes in a n x n (micrometers) neighborhood containing fewer than q bead barcodes
message("-------- Step 3: Remove bead barcodes in a ", n, " x ", n, " um neighborhood containing fewer than ", q, " bead barcodes")
embedding2<-seuratO@reductions$SPATIAL@cell.embeddings[grouping_filtered[,3],]
x <- cut(embedding2[,1], breaks = grid_n, labels = seq(1,grid_n))
y <- cut(embedding2[,2], breaks = grid_n, labels = seq(1,grid_n))
grouping <- data.frame(rownames(embedding2), x, y)

groups <- grouping[,c("x","y")]
#library(dplyr)
density <- dplyr::count(groups%>%group_by_all)

png(file.path(datadir, paste0(sample_id, "_step3_density_", n,"_", q,".png")), width=4.5, height=4.5, res = 300, units = "in")
plot(density(density$n), xlab=paste0("# bead barcodes / (",n,"x",n,") um"), main="")
dev.off()

density_high<-as.data.frame(density[density$n>=q,])
grouping_filtered <- merge(grouping, density_high, by=c("x", "y"))


message("-------- Visualize on/off tissue signals after step 3 filter")

## ON (#bead barcodes >= q on a nxn um neighborhood)
png(file.path(datadir, paste0(sample_id,"_step3_UMIgte", min_nUMI_m, "_ontissue_", m,"_",p,"_", n, "_", q, ".png")), width=7, height=4.5, res = 300, units = "in")
DimPlot(seuratO, reduction = "SPATIAL", cells.highlight = grouping_filtered[,3], cols= "white", raster=FALSE)&
  coord_fixed()&NoAxes()&
  theme(legend.position = "none")&
  ggtitle(paste0("Number of beads: " ,length(grouping_filtered[,3])))
dev.off()

## OFF (#bead barcodes < q on a nxn um neighborhood)
png(file.path(datadir, paste0(sample_id,"_step3_UMIgte", min_nUMI_m, "_offtissue_", m,"_",p,"_", n, "_", q, ".png")), width=7, height=4.5, res = 300, units = "in")
DimPlot(seuratO, reduction = "SPATIAL", cells.highlight = setdiff(grouping[,1], grouping_filtered[,3]), cols= "white", raster=FALSE)&
  coord_fixed()&NoAxes()&
  theme(legend.position = "none")&
  ggtitle(paste0("Number of beads: " ,length(setdiff(grouping[,1], grouping_filtered[,3]))))
dev.off()


# Output retained bead barcodes
message("-------- Step 3: Output retained bead barcodes")
write.table(grouping_filtered[,3], file.path(datadir, paste0(sample_id, "_UMIgte", min_nUMI_m, "_ontissue_", m,"_",p,"_", n, "_", q, ".txt")), row.names = FALSE, col.names = FALSE, quote = FALSE)

