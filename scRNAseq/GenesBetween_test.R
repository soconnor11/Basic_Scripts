

#------------------------
# Set up section / Load packages
#-----------------------

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
resdir = "update_classifier"

#------------------------------------------------------
# Load in NSC data
#---------------------------------------------------

tag = "hNSC"
data_dir <- 'ccAF/U5/filtered_gene_bc_matrices/hg19/'

# Genes as gene symbols
data <- Read10X(data.dir = data_dir, gene.column=2)
dim(data) #[1] 32738  3049

# Create seurat object
seurat1 <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
seurat1
dim(seurat1) #[1] 14734  3049

# Genes as gene symbols
mito.genes <- grep("MT-", rownames(seurat1))
percent.mito <- Matrix::colSums(seurat1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat1@assays[["RNA"]])
seurat1 <- AddMetaData(seurat1, percent.mito, col.name = "percent.mito")

pdf(file.path(resdir, paste0(tag, "_QC_plot.pdf")))
plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 2000, col = "red", lwd =3, lty =2)
abline(v = 130000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.1, col = "red", lwd =3, lty =2)
dev.off()

# Quality control filtering
keep.detect <- which(seurat1@meta.data$percent.mito < 0.1 & seurat1@meta.data$percent.mito > 0.009 &
                        seurat1@meta.data$nCount_RNA < 130000 & seurat1@meta.data$nCount_RNA > 2000)
length(keep.detect) #[1] 2980
ncol(seurat1) - length(keep.detect) # lose 69 cells
seurat1 <- subset(seurat1, cells=colnames(seurat1)[keep.detect])
dim(seurat1) #[1] 14734  2980

# Normalization
seurat1 <- SCTransform(seurat1, vars.to.regress = "nCount_RNA", verbose = FALSE)

# Regress out cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat1 <- CellCycleScoring(seurat1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat1 <- ScaleData(seurat1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat1))
dim(seurat1) #[1] 13870  2980

# Save as RDS object
saveRDS(seurat1 , file.path(resdir, paste0(tag, "_cc_regressed.rds")))


#------------------------------------------------------
# Load in hBMMSC data
#---------------------------------------------------

tag = "hBMMSC"
data_dir <- 'ccAF_newData/hBMMSC/filtered_feature_bc_matrix/'

# Genes as gene symbols
data <- Read10X(data.dir = data_dir, gene.column=2)
dim(data) #[1] 36601  7207

# Create seurat object
seurat2 <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
seurat2
dim(seurat2) #[1] 20421  7145

# Genes as gene symbols
mito.genes <- grep("MT-", rownames(seurat2))
percent.mito <- Matrix::colSums(seurat2@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat2@assays[["RNA"]])
seurat2 <- AddMetaData(seurat2, percent.mito, col.name = "percent.mito")

pdf(file.path(resdir, paste0(tag, "_QC_plot.pdf")))
plot(seurat2@meta.data$nCount_RNA, seurat2@meta.data$percent.mito,
     xlab = "nCount_RNA", ylab = "percent.mito", pch = 20)
abline(v = 6000, col = "red", lwd =3, lty =2)
abline(v = 130000, col = "red", lwd =3, lty =2)
abline(h = 0.009 , col = "red", lwd =3, lty =2)
abline (h = 0.1, col = "red", lwd =3, lty =2)
dev.off()

# Quality control filtering
keep.detect <- which(seurat2@meta.data$percent.mito < 0.1 & seurat2@meta.data$percent.mito > 0.009 &
                        seurat2@meta.data$nCount_RNA < 130000 & seurat2@meta.data$nCount_RNA > 6000)
length(keep.detect) # 6171
ncol(seurat2) - length(keep.detect) # lose 974 cells
seurat2 <- subset(seurat2, cells=colnames(seurat2)[keep.detect])
dim(seurat2) #[1] 20421  6171

# Normalization
seurat2 <- SCTransform(seurat2, vars.to.regress = "nCount_RNA", verbose = FALSE)

# Regress out cell cycle
seurat2 <- CellCycleScoring(seurat2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat2 <- ScaleData(seurat2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat2))
dim(seurat2) #[1] 19216  6171

# Save as RDS object
saveRDS(seurat2 , file.path(resdir, paste0(tag, "_cc_regressed.rds")))


#------------------------------------------------------
# Combine data
#---------------------------------------------------

seurat1 <- readRDS("update_classifier/hNSC_cc_regressed.rds")
dim(seurat1) #[1] 13870  2980
seurat2 <- readRDS("update_classifier/hBMMSC_cc_regressed.rds")
dim(seurat2) #[1] 19216  6171

data.combined <- merge(seurat1, y = seurat2, add.cell.ids = c("NSC", "BMMSC"))
data.combined

#DefaultAssay(object = data.combined) <- "RNA"
#data.combined_2 <- FindVariableFeatures(data.combined, selection.method = "vst")

seurat1 <- FindVariableFeatures(seurat1, selection.method = "vst")
seurat2 <- FindVariableFeatures(seurat2, selection.method = "vst")
lFeatCmpr <- list(
  S1 = seurat1@assays$RNA@meta.features,
  S2 = seurat2@assays$RNA@meta.features
)
lFeatCmpr <- lapply(lFeatCmpr, function(a) {
  a$Genes <- rownames(a)
  return(a)
})
lFeatCmprDf <- reshape2::melt(lFeatCmpr, id.vars = colnames(lFeatCmpr[[1]]))
meanCmpr <- lFeatCmprDf %>% reshape2::dcast(Genes~L1, value.var = "vst.mean")
varCmpr <- lFeatCmprDf %>% reshape2::dcast(Genes~L1, value.var = "vst.variance")

threshDf <- mutate(meanCmpr, AbsDiff = abs(S1 - S2))
DiffGenesBetween <- subset(na.omit(threshDf), AbsDiff > quantile(AbsDiff, 0.75)) %>% pull(Genes)
lFeatCmprDf2 <- merge(lFeatCmprDf, dplyr::select(na.omit(threshDf), Genes, AbsDiff))

HighBetweenLowWithinDf <- subset(lFeatCmprDf2, vst.variance < AbsDiff/4)
HighBetweenLowWithinGenes <- HighBetweenLowWithinDf %>% count(Genes) %>% subset(n == 2) %>% pull(Genes)

SubFeatCmprDf <- subset(lFeatCmprDf, Genes %in% DiffGenesBetween)
HighBetweenLowWithinGenes <- subset(SubFeatCmprDf, vst.variance < 0.1) %>% count(Genes) %>% subset(n == 2) %>% pull(Genes) # modify threshold as you need
