
################################################################
## Cell Deconvolution with Nowakowski_2017 - Data Preparation ##
################################################################

# Samantha O'Connor
# School of Biological Health Systems and Engineering
# Plaisier Lab
# Arizona State University
# saoconn1@asu.edu
# Last updated: September 14, 2022

#############################################

docker run -it -v '/home/soconnor/old_home/scGlioma/Nowakowski/:/files' saoconn1/scrnaseq_latest

#--------------------------------
# Set up section / load packages
#--------------------------------

# Load packages
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
library(reshape2)

# Set working directory
setwd("files/")

# Create folder to save out figures / csv files / etc.
tag = 'Nowakowski'
dirs = c('analysis_output', 'final_output')
for (dir in dirs){
  dir.create(dir, showWarnings = FALSE)
}
resdir = dirs[1]
resdir1 = dirs[2]


#------------------------------------------------------
# Read in Nowakowski data / set up seurat object
#---------------------------------------------------

seurat1 <- readRDS(file.path(resdir, paste0(tag, "_unfiltered.rds")))

#------------------------------------------------------
# Remove NAs & filter to cells from all areas of brain except for "0"
#---------------------------------------------------

# Remove NAs
RowsNA <- names(which(is.na(seurat1$WGCNAcluster)))
cells <- colnames(seurat1)[!(colnames(seurat1) %in% RowsNA)]
subset1 <- subset(seurat1, cells = cells)
sum(is.na(subset1@meta.data))
dim(subset1) #[1] 41222  4129

# Remove "0" Area
aoi = c('PFC', 'M1', 'V1', 'MGE', 'LGE', 'Choroid Plexus')
Idents(object = subset1) <- "Area"
subset1 <- subset(subset1, idents = aoi)
dim(subset1) #[1] 41222  3047

# Save out filtered RDS file
saveRDS(subset1, file.path(resdir1, paste0(tag, "_filtered.rds")))

#------------------------------------------------------
# Downstream processing
#---------------------------------------------------

# Load in filtered seurat object
#subset1 <- readRDS(file.path(resdir2, paste0(tag, "_filtered.rds")))

# save as new object so can go back to previous non-normalized / scaled seurat object if need to
seurat2 <- subset1
seurat2 <- SCTransform(seurat2, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
seurat2 <- RunPCA(seurat2, dims = 1:30)
seurat2 <- FindNeighbors(seurat2, dims = 1:30)
seurat2 <- FindClusters(seurat2)
seurat2 <- RunUMAP(seurat2, dims=1:30)

# Save out normalized RDS file
saveRDS(seurat2, file.path(resdir1, paste0(tag, "_normalized.rds")))

# Visualize clusters; adjust above parameters if need be
m1 = DimPlot(seurat2, reduction = "umap", label=T, label.size = 6)
m2 = DimPlot(seurat2, reduction = "umap", group.by = 'Area')
m3 = DimPlot(seurat2, reduction = "umap", group.by = 'WGCNAcluster')

# Plot
pdf(file.path(resdir1, paste0(tag, "_umap_test.pdf")), height = 8, width = 10)
m1
m2
m3
dev.off()


#------------------------------------------------------
# Combine all cell types to cell type groups
#---------------------------------------------------

Idents(seurat2) <- seurat2$WGCNAcluster

# Group cell type by main cell groups
seurat3 <- RenameIdents(seurat2, "tRG" = "Glial_Progenitor", "IPC-nEN1" = "Excitatory", "IPC-div2" =  "Neuronal_Progenitor","IPC-div1" = "Neuronal_Progenitor", "nEN-early2" = "Excitatory", "OPC" = "Glial_Progenitor", "IN-CTX-MGE2" = "Inhibitory", "nEN-late" = "Excitatory", "Mural" = "Mural", "IPC-nEN3" = "Excitatory", "RG-div1" = "Glial_Progenitor","IN-CTX-MGE1" = "Inhibitory", "IN-CTX-CGE2" = "Inhibitory","oRG" = "Glial_Progenitor","RG-div2" = "Glial_Progenitor","vRG" = "Glial_Progenitor", "EN-V1-2" =  "Excitatory", "EN-V1-3" =  "Excitatory", "IN-CTX-CGE1" = "Inhibitory", "IN-STR" = "Inhibitory","IPC-nEN2" = "Excitatory","EN-PFC3" = "Excitatory", "MGE-RG2" = "Glial_Progenitor", "nIN5" = "Inhibitory", "MGE-IPC3" = "Neuronal_Progenitor", "nIN4" = "Inhibitory", "MGE-IPC1" = "Neuronal_Progenitor", "MGE-RG1" = "Glial_Progenitor", "MGE-div" = "Neuronal_Progenitor", "MGE-IPC2" = "Neuronal_Progenitor", "nIN2" = "Inhibitory", "nIN1" = "Inhibitory","EN-V1-1" =  "Excitatory", "nIN3" =  "Inhibitory", "EN-PFC1" =  "Excitatory", "nEN-early1" =  "Excitatory", "EN-PFC2" =  "Excitatory", "RG-early" = "Glial_Progenitor", "Endothelial" = "Endothelial", "Astrocyte" = "Astrocyte", "Microglia" = "Microglia", "Glyc" = "Glyc",  "U1" = "U1", "Choroid" = "Choroid")

# Add new identities to a column in the object & save out new identities as a csv (because idents don't save out when you read in the RDS file, so you have to add back into the object if reading in)
seurat3$cell_type_combined <- Idents(seurat3)
write.csv(data.frame(seurat3$cell_type_combined), file.path(resdir1, paste0(tag, "_combined_cell_types.csv")))


#------------------------------------------------------------------
# Save as loom file for prepare_for_Cibersortx.py script (python)
#------------------------------------------------------------------

# Read back in to save out as loom file / proceed with script if starting from here
seurat3 <- readRDS(file.path(resdir1, paste0(tag, "_normalized.rds")))
tmp <- read.csv(file.path(resdir1, paste0(tag, "_combined_cell_types.csv")))
rownames(tmp) <- tmp$X
seurat3$cell_type_combined = tmp$seurat3.cell_type_combined


data_loom <- as.loom(seurat3, file.path(resdir1, paste0(tag, "_normalized_celltypes_combined.loom")), verbose = FALSE)
data_loom$close_all()


#-------------------------------------------------------------------------
# Find marker genes for each cell type - important for medioid generation / deconvolution
#-----------------------------------------------------------------------

Idents(seurat3) <- seurat3$cell_type_combined
cluster_markers= FindAllMarkers(seurat3, logfc.threshold = 0.25, only.pos = TRUE) # positive marker genes
cluster_no = unique(cluster_markers$cluster)
write.csv(cluster_markers, file.path(resdir1, paste0(tag,"_scTransform_Markers_together_celltypes_combined.csv")))
