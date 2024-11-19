setwd('/Users/cplaisie/Dropbox (ASU)/TBI Sex hormone data shared/')

# Load counts
d0 = read.csv('data_gene_count.csv',header=T,row.names=1)
d0 = d0[,1:33]
d0 = d0[rowSums(d0)!=0,]

library(DESeq2)

# Structure data
#conds = data.frame(id=colnames(d0))
conds = t(read.csv('meta_data.csv',header=T,row.names=1))

# Write out normalized data
dds1 = DESeqDataSetFromMatrix(countData = d0, colData = conds, design= ~ 1) # subtype
rld1 = rlog(dds1,blind=T)
write.csv(assay(rld1),'gexp_norm.csv')


##################
#### Analyses ####
##################

## Clustering samples
library(pheatmap)
pdf('correlationReplicates.pdf')
#c1 = cor(assay(rld)[names(sort(apply(assay(rld), 1, mad),decreasing=T))[1:10000],],method='spearman')
c1 = cor(assay(rld1),method='pearson')
diag(c1) = NA
pheatmap(c1,color=colorRampPalette(c('blue4','white','gold'))(32)) # ,clustering_method='single') # Looks much better
dev.off()


## PCA plots
# For groups with 12 or less samples
plotPCA_v2 <- function(x, intgroup = "condition", ntop = 500, col)
{
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]),
        1, paste, collapse = " : "))
    if (missing(col)) {
        col = if (nlevels(fac) >= 3)
            brewer.pal(nlevels(fac), "Paired")
        else c("red", "blue")
    }
    imp = summary(pca)$importance[2,]
    par(mar=c(1,1,1,1))
    xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x),
        pch = 16, cex = 2, aspect = "iso", col = col, main = draw.key(key = list(rect = list(col = col),
            text = list(levels(fac)), rep = FALSE)), xlab = paste('PC1 (',signif(imp[[1]],2)*100,'%)',sep=''), ylab = paste('PC2 (',signif(imp[[2]],2)*100,'%)',sep=''))
}

# Some libraries we need to run the PCA anlayses
library(matrixStats)
library(RColorBrewer)
library(lattice)

# Do the actual ploting
pdf('pca_rlog_counts.pdf')
plotPCA_v2(rld1, intgroup=c('Sex','Group'),ntop=2000)
plotPCA_v2(rld1, intgroup=c('Phase'),ntop=2000)
dev.off()
