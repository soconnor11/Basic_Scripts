
# Guided tutorial if need help / suggestions
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

#ssh soconnor@plaisierlab1.sbhse.dhcp.asu.edu
#docker run -it -v '/home/soconnor/ccAF_newData:/files' soconnor/seurat:scrnaseq

#--------------------------------
# Set up section / load packages
#--------------------------------

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
sessionInfo()
library(ssgsea.GBM.classification)
library(verification)
library(MCMCpack)


#------------------------------------------------------
# Read in data section / set up seurat object / QC
#---------------------------------------------------

# Set working directory
setwd("files/")
resdir = "redo_analysis"


#------------------------------------------------------
# Load in data (10X)
#---------------------------------------------------

tag = 'hBMMSC' # useful for saving figures
data_dir <- 'hBMMSC/filtered_feature_bc_matrix/' # where your data is located


#------ USE ENSEMBL IDs IF WANT TO SAVE OUT AS LOOM FILE TO APPLY CCAF --------#
# Genes as ensembl IDs - export as loom file to apply ccAF classification
#data <- Read10X(data.dir = data_dir, gene.column=1) # column 1 is ensembl IDs (in 10X mtx file)
#dim(data) #[1] 36601  7207
#------------------------------------------------------------------------------#

# Load in data
# Genes as gene symbols
data <- Read10X(data.dir = data_dir, gene.column=2) # column 2 is gene symbols (in 10X mtx file)
dim(data)

# Create seurat object
seurat1 <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
seurat1
dim(seurat1)


#------------------------------------------------------
# Quality control
#---------------------------------------------------

#------ USE ENSEMBL IDs IF WANT TO SAVE OUT AS LOOM FILE TO APPLY CCAF --------#
#mitogenes <- read_csv("hBMMSC/mito_genes.csv")
#list = as.list(mitogenes)
#seurat1[["percent.mito"]] <- PercentageFeatureSet(seurat1, features = as.list(mitogenes)$mito)/100


#------ if have ENSEMBL IDs  -------------------#
mitogenes <- read_csv("mito_genes.csv")
list = c("ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712", "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000198840", "ENSG00000212907", "ENSG00000198886", "ENSG00000198786", "ENSG00000198695", "ENSG00000198727")

list2 = c("ENSG00000210049", "ENSG00000211459", "ENSG00000210077", "ENSG00000210082", "ENSG00000209082", "ENSG00000198888", "ENSG00000210100", "ENSG00000210107", "ENSG00000210112", "ENSG00000198763", "ENSG00000210117", "ENSG00000210127", "ENSG00000210135", "ENSG00000210140", "ENSG00000210144", "ENSG00000198804", "ENSG00000210151", "ENSG00000210154", "ENSG00000198712", "ENSG00000210156", "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000210164", "ENSG00000198840", "ENSG00000210174", "ENSG00000212907", "ENSG00000198886", "ENSG00000210176", "ENSG00000210184", "ENSG00000210191", "ENSG00000198786", "ENSG00000198695", "ENSG00000210194", "ENSG00000198727", "ENSG00000210195", "ENSG00000210196")
seurat1[["percent.mito"]] <- PercentageFeatureSet(seurat1, features = list)
list %in% rownames(seurat1)
list2 %in% rownames(seurat1)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


# Genes as gene symbols
# Find mitochondiral genes
mito.genes <- grep("MT-", rownames(seurat1))
percent.mito <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")

# Plot & adjust cutoffs before filtering seurat object
pdf(file.path(resdir, paste0(tag, "_QC_plot.pdf")))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 6000, col = "red", lwd =3, lty =2)
abline(v = 130000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.1, col = "red", lwd =3, lty =2)
dev.off()

# Before filtering data - look at above plot and adjust cutoffs!
# Once find optimal cutoffs, enter values into keep.detect

# Quality control filtering
keep.detect <- which(seurat1@meta.data$percent.mito < 0.1 & seurat1@meta.data$percent.mito > 0.009 &
                        seurat1@meta.data$nCount_RNA < 130000 & seurat1@meta.data$nCount_RNA > 6000)
length(keep.detect) # will give you amount of cells after filtering
ncol(seurat1) - length(keep.detect) # will tell you how many cells you filtered out with above cutoffs
seurat1 <- subset(seurat1, cells=colnames(seurat1)[keep.detect]) # actually filter
dim(seurat1) # gives you dimensions of filtered object

# Check if filtering worked - cells should be inside red line cutoffs (if outside, you have your < > backwards!)
# Make sure to match cutoff values with filtered values
pdf(file.path(resdir, paste0(tag, "_QC_plot_POST.pdf")))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 6000, col = "red", lwd =3, lty =2)
abline(v = 130000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.1, col = "red", lwd =3, lty =2)
dev.off()


#------------------------------------------------------
# Save out as loom file (check to make sure genes are ensembl IDs)
# Comment out when finished
#---------------------------------------------------

#data_loom <- as.loom(seurat1, filename = paste0(tag, "_data.loom"), verbose = FALSE)
#data_loom
#data_loom$close_all()


#------------------------------------------------------
# Downstream processing
#---------------------------------------------------

# Load ccAF classifications and add as column into data
#ccAF <- read_csv("hBMMSC_data_ccAF_calls.csv")
#seurat1 <- AddMetaData(object = seurat1, metadata = ccAF$ccAF, col.name = 'ccAF')


#------------------------------------------------------
# Clustering / seurat pipeline with fastMNN -scTransform version.
#---------------------------------------------------

# save as new object so can go back to previous non-normalized / scaled seurat object if need too
seurat2 <- seurat1

seurat2 <- SCTransform(seurat2, vars.to.regress = "nCount_RNA", verbose = FALSE)
seurat2 <- RunPCA(seurat2, dims = 1:30)
#seurat2 <- RunTSNE(seurat2, dims = 1:30)
seurat2 <- RunUMAP(seurat2,   dims = 1:30)
seurat2 <- FindNeighbors(seurat2,   dims = 1:30)
seurat2 <- FindClusters(seurat2)

#res = 0.8 # Choose resolution
#seurat2 <- FindClusters(seurat2, resolution = res)

m1 = DimPlot(seurat2, reduction = "umap", label=T, label.size = 6) # if you use TSNE, replace "umap" with "tsne" in reduction parameter

# Plot
pdf(file.path(resdir, paste0(tag, "_umap_test.pdf")), width = 8, height = 8)
m1
dev.off()

# Open figure just plotted & look to see if number of clusters makes sense with your data!! Literature search / look at original paper / etc. Looking at marker genes to find meaning for clusters also helps!
# If not, un-comment line 147 and change resolution (lower value, less clusters; higher value, more clusters). Then uncomment 148 and run.

#---------------------------------------------------
# Find marker genes and save to excel file
#---------------------------------------------------

cluster_markers= FindAllMarkers(seurat2, logfc.threshold = 0.25, only.pos = TRUE) # positive marker gens
table(cluster_markers$cluster)
cluster_no = unique(cluster_markers$cluster)

sp = split(cluster_markers, cluster_markers$cluster)
names(sp) = paste0("cluster-", names(sp))
write_xlsx(sp, path = file.path(resdir, paste0(tag, "_scTransform_Markers_by_cluster.xlsx")))

cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

print(top10)

#-------------------------------#
# When decide on a resolution:
#-------------------------------#

# Visualize UMAP
m1 = DimPlot(seurat2, reduction = "umap", label=T, label.size = 6) + ggtitle ("de novo clustering")
#af = DimPlot(seurat2, reduction = "umap", group.by="ccAF", label=T) + ggtitle("ccAF")


#------------------------------------------------------
# Make Violin Plots for each cluster
#---------------------------------------------------

cluster_no = length(unique(seurat2$seurat_clusters))
o1 = VlnPlot(object = seurat2, features= "nFeature_RNA",   pt.size=0.1 )
o2 = VlnPlot(object = seurat2, features= "nCount_RNA",  pt.size=0.1 )
o4 = VlnPlot(object = seurat2, features= "percent.mito",   pt.size=0.1)


#---------------------------------------------------
# cell cycle scoring - ccSeurat
#---------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat2 <- CellCycleScoring(object=seurat2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)

c1 =  DimPlot(seurat2, reduction = "umap", label=T, group.by="Phase") + ggtitle("ccSeurat")


#---------------------------------------------------
# save as RDS object and write plots to pdf
#---------------------------------------------------
saveRDS(seurat2 , file.path(resdir, paste0(tag, "_seurat2.rds")))

# Dump out pdf of plots
pdf(file.path(resdir, paste0(tag, "_test.pdf")), width =8, height = 10)

# umaps ( with ccAF cell cycle scoring )
#lst = list(m1, c1, af )
#grid.arrange(grobs = lst, layout_matrix = rbind(c(1,NA), c(2,NA), c(3,NA)), top = "")

lst = list(m1, c1 ) # no ccAF
grid.arrange(grobs = lst, layout_matrix = rbind(c(1,NA), c(2,NA)), top = "")

lst2 = list( o1, o2, o4)
grid.arrange(grobs = lst2, layout_matrix = rbind(c(1,2), c(3)), top = "")

# Heatmap
DoHeatmap(object = seurat2, features = top10$gene)

# Cylin ridgeplots
features <- c("CCND1", "CCNE2", "CCNA2", "CCNB1","CDK1","CDK2")
RidgePlot(seurat2, features, ncol=2)

# Stacked barplot
barplot(table(seurat2$Phase, seurat2$seurat_clusters), beside = FALSE, col = c("blue", "green", "red"), xlab = "clusters ID", ylab = "Cell count", legend.text = rownames(table(seurat2$Phase, seurat2$seurat_clusters)), args.legend=list(title="ccSeurat classification"))

dev.off()


#---------------------------------------------------
# Assigning cell type identity to clusters
#---------------------------------------------------

#new.cluster.ids <- c("celltype1", "celltype2", "celltype3") # number of new cluster IDs needs to be same as number of leiden clusters (also in same order - leiden cluster 0 = "celltype1" in this list and so on...)
new.cluster.ids <- c("Late G1", "M", "M/Early G1", "G1", "S/G2", "G2",
    "G1/other", "G0", "S")
names(new.cluster.ids) <- levels(seurat2)
seurat2 <- RenameIdents(seurat2, new.cluster.ids) # rename clusters in seurat object
seurat2$new_clusters <- seurat2@active.ident
m2 = DimPlot(seurat2, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4) + ggtitle("labeled clusters")

# Plot with labeled umap
pdf(file.path(resdir, paste0(tag, "_labeled_test.pdf")), width =8, height = 8)

# umaps
lst3 = list(m1, m2 )
grid.arrange(grobs = lst3, layout_matrix = rbind(c(1,NA), c(2,NA)), top = "")

# Stacked barplot
barplot(table(seurat2$Phase, seurat2$new_clusters), beside = FALSE, col = c("blue", "green", "red"), xlab = "new cluster classification", ylab = "Cell count", legend.text = rownames(table(seurat2$Phase, seurat2$new_clusters)), args.legend=list(title="ccSeurat classification"))

dev.off()


#---------------------------------------------------
# Additional resources!!!
#---------------------------------------------------

# use Seurat API and PBMC tutorial for help with paramater functions & plot examples!
# https://satijalab.org/seurat/
