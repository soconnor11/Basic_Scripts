# docker run -it -v '/home/soconnor/scGlioma:/files' cplaisier/scrna_seq_velocity
​# cd /files
#R

# adapted functions (new parameters) in CONICSmat
#cd CONICS
#R
#devtools::install('./CONICSmat')
#1
#ctr d
#cd ..
#R

library(CONICSmat)

tag = "GSC827_inVitro_combined"
#suva_expr = as.matrix(read.table("OG_processed_data_portal.txt",sep="\t",header=T,row.names=1,check.names=F))
suva_expr = as.matrix(read.table(paste0(tag,"_OG_processed_data_portal.txt") ,sep="\t",header=T,row.names=1,check.names=F))
suva_expr [which(is.na(suva_expr ))]=0
dim(suva_expr)
suva_expr[1:5,1:5]


# Load cluster information
#clusters = read.csv(paste0(tag,"_cluster_identities.csv"))
clusters = read.csv(paste0(tag,"_cluster_identities_seurat_clusts.csv"))
head(clusters)
#clusters1 = clusters$ccAF
clusters1 = clusters$seurat_clusters
clusters1 = as.character(clusters1)

#patients=unlist(strsplit(colnames(suva_expr),"_",fixed=TRUE))[seq(1,(3*ncol(suva_expr))-1,3)]
#unique(patients)
#patients[which(patients=="93")]="MGH93"
#patients[which(patients=="97")]="MGH97"

regions=read.table("chromosome_arm_positions_grch38.txt",sep="\t",row.names = 1,header = T)
head(regions,n=5)

# Preprocessing
gene_pos=getGenePositions(rownames(suva_expr))
#gene_pos=getGenePositions_ensembl(rownames(suva_expr)
#suva_expr=filterMatrix(suva_expr,gene_pos[,"hgnc_symbol"],minCells=5)
suva_expr=filterMatrix(suva_expr,gene_pos[,"hgnc_symbol"],minCells=5)
normFactor=calcNormFactors(suva_expr)
l=plotAll(suva_expr,normFactor,regions,gene_pos,tag)
# save out, read into scanpy, push into an obs and paint
write.csv(l, paste0(tag,"_posterior_probabilities_CONICSmap.csv"))
# python load into pandas as dataframe
# obs scanpy, axis=1

hi=plotHistogram(l,suva_expr,clusters=1,zscoreThreshold=4,clusters1)
# 20q specific to cells - paint that onto UMAPs
dev.off()


vg=detectVarGenes(suva_expr,500)

ts=calculateTsne(suva_expr,vg)
#plotTsneGene(ts,suva_expr,c("MBP","CSF1R","ALDOC","OLIG1"))

plotTsneProbabilities(ts,suva_expr,l[,"1p"],"Tsne Chr1p")

lrbic=read.table("SUVA_CNVs_BIC_LR.txt",sep="\t",header=T,row.names=1,check.names=F)
colnames(lrbic)
candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>200 & lrbic[,"LRT adj. p-val"]<0.01)]
hi=plotHistogram(l[,candRegions],suva_expr,clusters=4,zscoreThreshold=4,patients)
