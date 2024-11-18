
# Genes of interest
genes = c('ENSG00000106069','ENSG00000184672', 'ENSG00000185774')
# Set up dataframe for first orig.ident ('GF')
df_gf = data.frame(matrix(ncol=5, nrow=6))
colnames(df_gf) = c('cluster', genes, 'extra')
# Get expression values for genes of interest per identity
d1 = t(as.matrix(AverageExpression(seurat2, assays = 'RNA', features = genes, slot = "data", group.by = "seurat_clusters", verbose = FALSE)$RNA))
# Fill in dataframe
df_gf$cluster = rownames(d1)
for(gene1 in genes){
  df_gf[[gene1]] = d1[,gene1]
}
df_gf$extra = 'GF'

# Set up dataframe for second orig.ident ('noGF')
df_nogf = data.frame(matrix(ncol=5, nrow=6))
colnames(df_nogf) = c('cluster', genes, 'extra')
# Get expression values for genes of interest per identity
d2 = t(as.matrix(AverageExpression(seurat2, assays = 'RNA', features = genes, slot = "data", group.by = "seurat_clusters", verbose = FALSE)$RNA))
# Fill in dataframe
df_nogf$cluster = rownames(d2)
for(gene1 in genes){
  df_nogf[[gene1]] = d2[,gene1]
}
df_nogf$extra = 'noGF'

# Combine dataframes
df_combine = rbind(df_gf, df_nogf)
df_combine_reshape = reshape2::melt(df_combine, id = c("cluster","extra"))

# Plot
ggplot(df_combine_reshape, aes(x = extra, y = value, fill = variable)) + geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ cluster)
