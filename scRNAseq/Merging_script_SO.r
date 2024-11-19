
docker run -it -v '/home/soconnor:/files' soconnor/seurat:scrnaseq

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
library(stringr)


#----------------------------
# Set up working directory
#----------------------------

# Set working directory
setwd("files/rawls/")
resdir <- "redo_analysis"

#------------------------------------------------------
# Load in data
#---------------------------------------------------

dirs = c('WT_L', 'MKX_L')
tags = c('WT_L', 'MKX_L')
tag = tags[1]
dir = dirs[1]

# WT_L
wt <- Read10X(data.dir = dir, gene.column=2)
dim(wt) #[1] 32285  3021

# Create seurat object
seurat1 <- CreateSeuratObject(counts = wt, min.cells = 3, min.features = 200)
seurat1
dim(seurat1) #[1] 16557  2999

# adding mito meta data to each
mito.genes <- grep("mt-", rownames(seurat1))
percent.mito <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")

# ADJUST VERTICAL AND HORIZONTAL CUTOFF LINES DEPENDING ON DATA
v1 <- 2000
v2 <- 80000
h1 <- 0.009
h2 <- 0.2

pdf(file.path(resdir, paste0(tag, "_QC_plot.pdf")))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = v1, col = "red", lwd =3, lty =2)
abline(v = v2, col = "red", lwd =3, lty =2)
abline(h = h1 , col = "red", lwd =3, lty =2)
abline (h = h2, col = "red", lwd =3, lty =2)
dev.off()

# QC
keep.detect <- which(seurat1@meta.data$percent.mito < h2 & seurat1@meta.data$percent.mito > h1 &
                       seurat1@meta.data$nCount_RNA < v2 & seurat1@meta.data$nCount_RNA > v1)
length(keep.detect) #[1] 2751
ncol(seurat1) - length(keep.detect) # [1] 248
seurat1.1 <- subset(seurat1, cells=colnames(seurat1)[keep.detect])
dim(seurat1.1) #[1] 16557  2751
seurat1.1$type = tag

saveRDS(seurat1.1, file.path(resdir, paste0(tag, "_filtered.rds")))


#-------------------------------------------------------
  #transforming and normalizing

res = 0.45 # Choose resolution

seurat1.2 <- SCTransform(seurat1.1, vars.to.regress = "nCount_RNA", verbose = F)
seurat1.2 <- RunPCA(seurat1.2, verbose = FALSE)
seurat1.2 <- FindNeighbors(seurat1.2,   dims = 1:10, verbose = FALSE)
seurat1.2 <- FindClusters(seurat1.2, verbose = FALSE, resolution = res) # default resolution
seurat1.2 <- RunTSNE(seurat1.2)
seurat1.2 <- RunUMAP(seurat1.2, dims=1:10)

# test plots
m1 = DimPlot(seurat1.2, reduction = "umap", label=T, label.size = 6) + ggtitle ("de novo clustering umap")
m2 = DimPlot(seurat1.2, reduction = "tsne", label=T, label.size = 6) + ggtitle ("de novo clustering tsne")

pdf(file.path(resdir, paste0(tag, "_test.pdf")))
m1
m2
dev.off()

saveRDS(seurat1.2, file.path(resdir, paste0(tag, "_normalized.rds")))
write.csv(seurat1.2$seurat_clusters, file.path(resdir, paste0(tag, "_seurat_cluster_labels.csv")))



#---- MKX_L-----#

tag = tags[2]
dir = dirs[2]

mkx = Read10X(data.dir = dir, gene.column=2)
dim(mkx)
seurat2 = CreateSeuratObject(counts = mkx, min.cells = 3, min.features = 200)
dim(seurat2)

mito.genes <- grep("mt-", rownames(seurat2))
percent.mito <- Matrix::colSums(seurat2@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat2@assays[["RNA"]])
seurat2 <- AddMetaData(seurat2, percent.mito, col.name = "percent.mito")

# ADJUST VERTICAL AND HORIZONTAL CUTOFF LINES DEPENDING ON DATA
v1 <- 2000
v2 <- 80000
h1 <- 0.009
h2 <- 0.2

pdf(file.path(resdir, paste0(tag, "_QC_plot.pdf")))
plot(seurat2@meta.data$nCount_RNA, seurat2@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = v1, col = "red", lwd =3, lty =2)
abline(v = v2, col = "red", lwd =3, lty =2)
abline(h = h1 , col = "red", lwd =3, lty =2)
abline (h = h2, col = "red", lwd =3, lty =2)
dev.off()


keep.detect2 <- which(seurat2@meta.data$percent.mito < h2 & seurat2@meta.data$percent.mito > h1 &
                        seurat2@meta.data$nCount_RNA < v2 & seurat2@meta.data$nCount_RNA > v1  )
length(keep.detect2) #[1] 1343
ncol(seurat2) - length(keep.detect2) # [1] 1877
seurat2.1 <- subset(seurat2, cells=colnames(seurat2)[keep.detect])
dim(seurat2.1) #[1] 16557  2688
seurat2.1$type = tag

saveRDS(seurat2.1, file.path(resdir, paste0(tag, "_filtered.rds")))


#-------------------------------------------------------
  #transforming and normalizing

res = 0.45 # Choose resolution

seurat2.2 <- SCTransform(seurat2.1, vars.to.regress = "nCount_RNA", verbose = F)
seurat2.2 <- RunPCA(seurat2.2, verbose = FALSE)
seurat2.2 <- FindNeighbors(seurat2.2,   dims = 1:10, verbose = FALSE)
seurat2.2 <- FindClusters(seurat2.2, verbose = FALSE, resolution = res) # default resolution
seurat2.2 <- RunTSNE(seurat2.2)
seurat2.2 <- RunUMAP(seurat2.2, dims=1:10)

# test plots
m1 = DimPlot(seurat2.2, reduction = "umap", label=T, label.size = 6) + ggtitle ("de novo clustering umap")
m2 = DimPlot(seurat2.2, reduction = "tsne", label=T, label.size = 6) + ggtitle ("de novo clustering tsne")

pdf(file.path(resdir, paste0(tag, "_test.pdf")))
m1
m2
dev.off()

saveRDS(seurat2.2, file.path(resdir, paste0(tag, "_normalized.rds")))
write.csv(seurat2.2$seurat_clusters, file.path(resdir, paste0(tag, "_seurat_cluster_labels.csv")))


#---------------------------------------------------------------------------------
# Set up working directory / Organize files to load - generate integrated dataset
#---------------------------------------------------------------------------------

# Set working directory
setwd("files/rawls/redo_analysis")

resdir <- "filtered_data"
files <- list.files(path = resdir, pattern = ".rds")

bm.list = c()
for (f1 in files){
  d1 = readRDS(file.path(resdir, f1))
  tmp <- f1
  tag <- str_split(tmp, "_")[[1]][1]
  d1$sample = tag
  lab = read.csv(paste0(tag, "_seurat_cluster_labels.csv"))
  d1 <- AddMetaData(object = d1, metadata = lab$x, col.name = 'prior_seurat_clusters')
  bm.list = append(bm.list,d1)
}

# normalize and identify variable features for each dataset independently
bm.list <- lapply(X = bm.list, FUN = function(x) {
    x <- SCTransform(x)
})


# Select features for downstream integration
bm.features <- SelectIntegrationFeatures(object.list = bm.list)
options(future.globals.maxSize = 8000 * 1024^5)
bm.list <- PrepSCTIntegration(object.list = bm.list, anchor.features = bm.features,
    verbose = FALSE)

# Identify anchors and integrate the datasets
bm.anchors <- FindIntegrationAnchors(object.list = bm.list, normalization.method = "SCT", anchor.features = bm.features, verbose = FALSE)
bm.integrated <- IntegrateData(anchorset = bm.anchors, normalization.method = "SCT", verbose = FALSE)

# Proceed with downstream analysis
bm.integrated_2 <- RunPCA(bm.integrated, verbose = FALSE)
bm.integrated_2 <- FindNeighbors(bm.integrated_2,   dims = 1:10, verbose = FALSE)
bm.integrated_2 <- FindClusters(bm.integrated_2, verbose = FALSE) # default resolution
bm.integrated_2 <- RunUMAP(bm.integrated_2, dims = 1:30)
plot1 = DimPlot(bm.integrated_2, combine = FALSE, label=T)
plot2 = DimPlot(bm.integrated_2, group.by = 'sample')


tag = 'integrated'
resdir2 = 'integrated'
# Plot

pdf(file.path(resdir2, paste0(tag, "_test.pdf")))
plot1
plot2
dev.off()

DefaultAssay(bm.integrated_2) <- "SCT" # important for looking at gene expression differences
plot3 = VlnPlot(bm.integrated_2, features= "Ccl2",   pt.size=0.1)

pdf(file.path(resdir2, paste0(tag, "_expression_test.pdf")))
plot3
dev.off()


# https://satijalab.org/seurat/articles/integration_introduction.html



#---------------------------------------------------
# identifing conserved cell type markers
#---------------------------------------------------
log2_thres = log2(1.25)
fdr_thres=0.05
cluster_markers = FindAllMarkers(merged1, logfc.threshold = 0.25, only.pos = TRUE) # positive marker gens
table(cluster_markers$cluster)
cluster_no = unique(cluster_markers$cluster)

sp = split(cluster_markers, cluster_markers$cluster)
names(sp) = paste0("cluster-", names(sp))
write_xlsx(sp, path = file.path(resdir, paste0(tag, "_scTransform_Markers_by_cluster_combined.xlsx")))

#ifn_conserved.markers = FindConservedMarkers(merged1, ident.1 = '12', grouping.var = "type", only.pos = T)
#head(ifn_conserved.markers)
#tail(ifn_conserved.markers)
#ifn_conserved.markers$gene = rownames(ifn_conserved.markers)
#write_xlsx(ifn_conserved.markers, path = file.path(resdir, paste0(tag, "_scTransform_Markers_by_cluster_conserved_markers.xlsx")))




#---------------------------------------------------
# Plotting dif ex
#---------------------------------------------------

markers.to.plot = c('Myod1','Pax7','Ccl2','Cxcl1','Cxcl10','Ccl5', 'Ccl7','Myog','Mef2c', 'Myf5', 'Myf6', 'Oasl1', 'Oasl2', 'Isg15')

d1 = DotPlot(merged1, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "type") + RotatedAxis()

d1

f1 = FeaturePlot(merged1, features =c('Myod1','Pax7','Ccl2','Cxcl1')  , col= c('gray', 'red'), split.by = 'type', keep.scale = 'feature')
f1
f2 = FeaturePlot(merged1, features =c('Cxcl10','Ccl5', 'Ccl7', 'Oasl1')  , col= c('gray', 'red'), split.by = 'type', keep.scale = 'feature')
f2
f3 = FeaturePlot(merged1, features =c('Myog','Mef2c', 'Myf5', 'Myf6','Myod1' )  , col= c('gray', 'red'), split.by = 'type', keep.scale = 'feature')
f3

#---------------------------------------------------
# Saving RDS object and outputting graphs
#---------------------------------------------------
saveRDS(merged1 , file.path(resdir, paste0(tag, "_data.rds")))



pdf(file.path(resdir, paste0(tag, "_output.pdf")), width =8, height = 10)
p1+p2
p3
v2
d1
f1
f2
f3

dev.off()
