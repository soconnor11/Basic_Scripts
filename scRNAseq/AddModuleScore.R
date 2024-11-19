
# Load in normalized object
seurat2 <- readRDS(file.path(resdir3, paste0(tag, "_normalized.rds")))
# Make sure the features are gene symbols
seurat2@assays
# Does this seurat object have dimensional reductions (umap) calculated? If not: Run downstream analysis
seurat2 <- RunPCA(seurat2, dims = 1:15)
seurat2 <- FindNeighbors(seurat2, dims = 1:15)
seurat2 <- FindClusters(seurat2, res = 0.5)
seurat2 <- RunUMAP(seurat2, dims=1:15)

# ccSeurat
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
#seurat2 <- CellCycleScoring(object=seurat2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)

# Read in integrated cluster information from paper (we changed to csv and deleted first row)
cellInfo = read.csv(file.path(resdir, '41467_2020_18075_MOESM2_ESM.csv'))
# Change all dashes to underscores (R doesn't like dashes)
cellInfo$Cluster <- sub("-", "_", cellInfo$Cluster)
# Remove rows where values were saved as dates in excel
cellInfo2 <- cellInfo[!cellInfo$Row %in% c("9-Sep", "10-Sep"),]

for(cluster1 in unique(cellInfo2$Cluster)){
  print(cluster1)
  # Find marker genes for specific clusters
  tmp = cellInfo2[cellInfo2$Cluster == cluster1,]
  tmp2 = list(tmp$Row)
  # Add gene modules to seurat object
  seurat2 = AddModuleScore(object = seurat2, features = tmp2, ctrl = 5, name = cluster1)
  # Remove "1" from column names so they are actually the cluster names from excel file
  seurat2[[cluster1]] = seurat2[[paste0(cluster1, '1')]]
  # Delete old column
  seurat2[[paste0(cluster1, '1')]] <- NULL
}

# Dim plots
d1 = DimPlot(seurat2, group.by = 'seurat_clusters', reduction = 'umap', size = 1) + ggtitle('de novo clusters')
d2 = DimPlot(seurat2, group.by = 'Phase', reduction = 'umap', size = 1) + ggtitle('ccSeurat classification')
d3 = DimPlot(seurat2, group.by = 'ccAF', reduction = 'umap', size = 1) + ggtitle('ccAF classification')
d4 = DimPlot(seurat2, group.by = 'ccAFv2', reduction = 'umap', size = 1) + ggtitle('ccAFv2 classification')

# Feature plots
f1 = FeaturePlot(seurat2, features = 'BAS_III', reduction = 'umap')
f2 = FeaturePlot(seurat2, features = 'MEL', reduction = 'umap')
f3 = FeaturePlot(seurat2, features = 'GRN', reduction = 'umap')
f4 = FeaturePlot(seurat2, features = 'BAS_IV', reduction = 'umap')
f5 = FeaturePlot(seurat2, features = 'BAS_II', reduction = 'umap')
f6 = FeaturePlot(seurat2, features = 'SPN', reduction = 'umap')
f7 = FeaturePlot(seurat2, features = 'BAS_I', reduction = 'umap')

# Plot all on one pdf
pdf(file.path(resdir2, 'integrative_cluster_from_paper.pdf'))
# Dim Plots
lst = list(d1, d2, d3, d4)
grid.arrange(grobs = lst, layout_matrix = rbind(c(1,2), c(3,4)), top = "")
# Feature Plots
lst1 = list(f1, f2, f3, f4)
grid.arrange(grobs = lst1, layout_matrix = rbind(c(1,2), c(3,4)), top = "")
lst2 = list(f5, f6, f7)
grid.arrange(grobs = lst2, layout_matrix = rbind(c(1,2), c(3,NA)), top = "")
# Violin Plot
VlnPlot(seurat2, features= c('BAS_III', 'MEL', 'GRN', 'BAS_IV', 'BAS_II', 'SPN', 'BAS_I'))
dev.off()
