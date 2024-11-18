#--------------------------------
# Set up section / load packages
#--------------------------------

library(dplyr)
library(Seurat)
library(SeuratData)
#library(SeuratDisk)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(readr)
sessionInfo()
#library(ssgsea.GBM.classification)
library(verification)
library(MCMCpack)
library(metap)
library(cowplot)

#------------------------------------------------------
# Read in data section / set up seurat object / QC
#---------------------------------------------------

# Set working directory
setwd("C:/Users/Alex/Dropbox (ASU)/scRNA-SEQ LPS wt mkx/4285926_Andre")
set
resdir = "analyzed looser gate med res CC regressed/WT_MKX_L/" 


#------------------------------------------------------
# Load in data
#---------------------------------------------------
tag = "WT_MKX_L_combined"
tag1 = "WT_L"
tag2 = 'Mkx_L'
data_dir.1 <- 'WT_L/outs/filtered_feature_bc_matrix/'
data.1 <- Read10X(data.dir = data_dir.1, gene.column=2)
dim(data.1) #[1] 32285  3021

# Create seurat object
seurat1 <- CreateSeuratObject(counts = data.1, min.cells = 3, min.features = 200)
seurat1
dim(seurat1) #[1] 16557  2999

data_dir.2 <- 'MKX_L/outs/filtered_feature_bc_matrix/'
data.2 = Read10X(data.dir = data_dir.2, gene.column=2)
dim(data.2) 
seurat2 = CreateSeuratObject(counts = data.2, min.cells = 3, min.features = 200)
dim(seurat2)

# adding mito meta data to each 


mito.genes.1 <- grep("mt-", rownames(seurat1))
percent.mito.1 <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes.1, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito.1, col.name = "percent.mito")



mito.genes.2 <- grep("mt-", rownames(seurat2))
percent.mito.2 <- Matrix::colSums(seurat2@assays[["RNA"]][mito.genes.2, ])/Matrix::colSums(seurat2@assays[["RNA"]])
seurat2 <- AddMetaData(seurat2, percent.mito.2, col.name = "percent.mito")



plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 2000, col = "red", lwd =3, lty =2)
abline(v = 80000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.2, col = "red", lwd =3, lty =2)

# QC
keep.detect <- which(seurat1@meta.data$percent.mito < 0.2 & seurat1@meta.data$percent.mito > 0.009 &
                       seurat1@meta.data$nCount_RNA < 80000 & seurat1@meta.data$nCount_RNA > 2000)
length(keep.detect) # 2688
ncol(seurat1) - length(keep.detect) # lose 311 cells
seurat1.1 <- subset(seurat1, cells=colnames(seurat1)[keep.detect])
dim(seurat1.1) #[1] 16557  2688
seurat1.1$type = tag1

plot(seurat2@meta.data$nCount_RNA, seurat2@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 2000, col = "red", lwd =3, lty =2)
abline(v = 80000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.2, col = "red", lwd =3, lty =2)


keep.detect2 <- which(seurat2@meta.data$percent.mito < 0.25 & seurat2@meta.data$percent.mito > 0.05 &
                        seurat2@meta.data$nCount_RNA < 130000 & seurat2@meta.data$nCount_RNA > 6000  )
length(keep.detect2) # 1073
ncol(seurat2) - length(keep.detect2) # lose 975 cells
seurat2.1 <- subset(seurat2, cells=colnames(seurat2)[keep.detect])
dim(seurat2.1) #[1] 16557  2688
seurat2.1$type = tag2

#-------------------------------------------------------
  #transforming and normalizing

res = 0.45 # Choose resolution

seurat1.2 <- SCTransform(seurat1.1, vars.to.regress = "nCount_RNA", verbose = T)
seurat1.2 <- RunPCA(seurat1.2, dims = 1:30)
seurat1.2 <- RunUMAP(seurat1.2,   dims = 1:30)
seurat1.2 <- FindNeighbors(seurat1.2,   dims = 1:30)
seurat1.2 <- FindClusters(seurat1.2, resolution=res)


seurat2.2 <- SCTransform(seurat2.1, vars.to.regress = "nCount_RNA", verbose = T)
seurat2.2 <- RunPCA(seurat2.2, dims = 1:30)
seurat2.2 <- RunUMAP(seurat2.2,   dims = 1:30)
seurat2.2 <- FindNeighbors(seurat2.2,   dims = 1:30)
seurat2.2 <- FindClusters(seurat2.2, resolution=res)






#------------------------------------------------------
# Merge seurat objects - use filtered seurat objects and normalized objects
#---------------------------------------------------

# Subset filtered Seurat objects based on common genes
common_genes <- intersect(rownames(seurat1.1), rownames(seurat2.1))
seurat1_com <- seurat1.1[common_genes,]
seurat2_com <- seurat2.1[common_genes,]

# Add normalized seurat object meta information to filtered (non-normalized) seurat objects
seurat1_com@meta.data <- seurat1.2@meta.data
seurat2_com@meta.data <- seurat2.2@meta.data

# Merge datasets

merged <- merge(seurat1_com, seurat2_com, add.cell.ids = c(tag1, tag2), project = "SC_LPS")
dim(merged) #[1] 18941  1152



 

# Normalize data together
merged1 <- SCTransform(merged, vars.to.regress = "nCount_RNA", verbose = T)
merged1 <- RunPCA(merged1, dims = 1:30)
merged1 <- RunUMAP(merged1,   dims = 1:30)
merged1 <- FindNeighbors(merged1,   dims = 1:30)
merged1 <- FindClusters(merged1, resolution=res)
m1 = DimPlot(merged1, reduction = "umap", label=T, label.size = 6) + ggtitle ("de novo clustering")
m1



# Visualization
p1 <- DimPlot(merged1, reduction = "umap", group.by = "type")
p2 <- DimPlot(merged1, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
p3 = DimPlot(merged1, reduction = "umap", split.by = "type")
p3

features =c('Myod1','Pax7','Ccl2','Cxcl1')
v2 = VlnPlot(merged1, features = features ,  split.by = 'type', pt.size = 0)
v2



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
