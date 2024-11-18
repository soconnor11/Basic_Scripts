
library(Seurat)
library(loomR)


#resdir = "U5_hNSC_Neural_G0/data/U5_hNSC/WT"
#tag = "U5_hNSC"
resdir = "scGlioma_1/GSE131928"
tag = "GSE131928"

# Import U5-NSC single cell data loom file (WT)
# Connect to the loom file in read/write mode
#lfile <- connect(file.path(resdir, "U5_all_v3.loom"), mode = "r+")
#lfile <- connect(file.path(resdir, "GSE131928_10X.loom"), mode = "r+")
lfile <- connect(file.path(resdir, "GSE131928.loom"), mode = "r+")
lfile

# Write as seurat object
seurat1 = as.Seurat(lfile)

# Save as RDS file
saveRDS(seurat1, file.path(resdir, paste0(tag, "_seurat1_10X_SmartSeq.rds")))
resdir2 = "ccAF_improve/gsea"
saveRDS(seurat1, file.path(resdir2, paste0(tag, "_seurat1.rds")))

# Save gene expression as txt for gsea classification
write.table(seurat1@assays$RNA@data, file.path(resdir2, paste0(tag, "_expression_data.txt")), sep="\t")





# Import U5-NSC single cell data loom file (integrated)
resdir = "U5_hNSC_Neural_G0/data"
resdir2 = "ccAF_improve/gsea"
tag = "integrated"
# Connect to the loom file in read/write mode
lfile <- connect(file.path(resdir, "cellcycle_int_integrated.loom"), mode = "r+")
lfile
# Write as seurat object
seurat1 = as.Seurat(lfile)
write.table(seurat1@assays$integrated@data, file.path(resdir2, paste0(tag, "_expression_data.txt")), sep="\t")
