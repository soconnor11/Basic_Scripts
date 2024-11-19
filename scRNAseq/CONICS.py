# docker run -it -v '/home/soconnor/gbmTriculture:/files' cplaisier/scrna_seq_velocity
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
"""
# if have problem in R can type name of function and change definiton (gene names to ensembl ID)

getGenePositions_ensembl = function (gene_names, ensembl_version = "dec2016.archive.ensembl.org",
    species = "human")
{
    if (species == "human") {
        ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "hsapiens_gene_ensembl", host = ensembl_version)
        gene_positions <- biomaRt::getBM(attributes = c("ensembl_gene_id",
            "hgnc_symbol", "chromosome_name", "start_position",
            "end_position"), filters = "ensembl_gene_id", values = gene_names,
            mart = ensembl)
    }
    else {
        ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "mmusculus_gene_ensembl", host = ensembl_version)
        gene_positions <- biomaRt::getBM(attributes = c("ensembl_gene_id",
            "mgi_symbol", "chromosome_name", "start_position",
            "end_position"), filters = "mgi_symbol", values = gene_names,
            mart = ensembl)
    }
    gene_positions = gene_positions[!duplicated(gene_positions[,
        2]), ]
    gene_positions[which(gene_positions[, 3] == "X"), 3] = 23
    gene_positions[which(gene_positions[, 3] == "Y"), 3] = 24
    gene_positions[which(gene_positions[, 3] == "MT"), 3] = 0
    gene_positions[which(nchar(gene_positions[, 3]) > 2), 3] = 0
    gene_positions = gene_positions[order(as.numeric(gene_positions[,
        3]), decreasing = F), ]
    return(gene_positions)
}


plotChrEnrichment_ensembl = function (expmat, chr, normFactor, gene_positions, n = 1, groups1 = NULL,
    groups2 = NULL, start = NULL, end = NULL, k = 2, vis = T,
    postProb = 0.95, repetitions = 5)
{
    par(mfrow = c(2, 2))
    if (!is.null(groups1)) {
        cellcolor = rep("black", (length(groups1) + length(groups2)))
        cellcolor[groups2] = "red"
    }
    else {
        cellcolor = NULL
    }
    if (!is.null(start)) {
        chr_genes = gene_positions[which(gene_positions[, 3] ==
            chr & gene_positions[, 4] > start & gene_positions[,
            5] < end), 1]
    }
    else {
        chr_genes = gene_positions[which(gene_positions[, 3] ==
            chr), 1]
    }
    if (length(chr_genes) > 100) {
        chr_exp = scale(colMeans(expmat[intersect(chr_genes,
            row.names(expmat)), ]) - normFactor)
        bestlog = (-Inf)
        bestmix = NULL
        loglik = NULL
        for (i in 1:repetitions) {
            print(paste("Fitting GMM for chr", chr, " ", start,
                ":", end, " iteration ", i, sep = ""))
            mixmdl = tryCatch(mixtools::normalmixEM(chr_exp,
                k = k, maxit = 1000, maxrestarts = 10), error = function(e) {
                print(paste("EM algorithm did not converge for region",
                  chr, " ", start, " ", end))
                mixmdl = NULL
            })
            if (!is.null(mixmdl)) {
                if (mixmdl$loglik > bestlog) {
                  bestlog = mixmdl$loglik
                  bestmix = mixmdl
                }
            }
        }
        if (is.null(bestmix)) {
            hist(chr_exp, breaks = 50, main = paste("Chr: ",
                chr, ":", start, ":", end, "\n", "Unable to fit 2 component mixture model",
                sep = ""), lwd2 = 3, xlab2 = "Expression z-score")
            plot(runif(length(chr_exp), 0, 100), chr_exp, pch = 16,
                ylab = "Expression z-score", ylim = c(min(chr_exp),
                  (max(chr_exp) + 2)), xlab = "Cells")
            hist(chr_exp, breaks = 50, main = paste("Chr: ",
                chr, ":", start, ":", end, "\n", "Unable to fit 2 component mixture model",
                sep = ""), lwd2 = 3, xlab2 = "Expression z-score")
            plot(runif(length(chr_exp), 0, 100), chr_exp, pch = 16,
                ylab = "Expression z-score", ylim = c(min(chr_exp),
                  (max(chr_exp) + 2)), xlab = "Cells")
        }
        else {
            out1 = list(x = chr_exp, mu = mean(chr_exp), sigma = sd(chr_exp),
                lambda = 1, loglik = sum(dnorm(chr_exp, mean(chr_exp),
                  sd(chr_exp), log = TRUE)))
            bics = c(max(BIC.mix(out1), 1), max(BIC.mix(bestmix),
                1))
            lrt = round(likelihoodRatioTest(out1$loglik, bestmix$loglik,
                n), 6)
            bestmix$BIC = bics
            bestmix$lrt = lrt
            if (vis == T) {
                plot(bestmix, which = 2, breaks = 50, col1 = c("red",
                  "green"), main2 = paste("Chr: ", chr, ":",
                  start, ":", end, "\n", "Log likelihood ", round(bestmix$loglik,
                    1), sep = ""), lwd2 = 3, xlab2 = "Expression z-score")
            }
            if (length(cellcolor) > 1 & vis == T) {
                g1 = length(which(bestmix$posterior[groups1,
                  1] > postProb))/length(groups1) * 100
                g2 = length(which(bestmix$posterior[groups1,
                  2] > postProb))/length(groups1) * 100
                g3 = length(which(bestmix$posterior[groups1,
                  2] < postProb & bestmix$posterior[groups1,
                  1] < postProb))/length(groups1) * 100
                g4 = length(which(bestmix$posterior[groups2,
                  1] > postProb))/length(groups2) * 100
                g5 = length(which(bestmix$posterior[groups2,
                  2] > postProb))/length(groups2) * 100
                g6 = length(which(bestmix$posterior[groups2,
                  2] < postProb & bestmix$posterior[groups2,
                  1] < postProb))/(length(groups2) + length(groups1)) *
                  100
                barplot(rbind(c(g1, g2, g3), c(g4, g5, g6)),
                  ylim = c(0, 100), beside = T, ylab = "Percentage of cells",
                  names = c("Cluster", "Cluster", "Ambigu"),
                  legend = c("Non-malignant", "Malignant"), args.legend = list(title = "Pred. via transcript.",
                    x = "topright", cex = 0.65), xlab = "Predicted via transcriptomics")
                axis(1, at = c(0.5, 1, 2, 3, 3.3), line = 2,
                  tick = T, labels = rep("", 5), lwd = 3, lwd.ticks = 0,
                  col = "red")
                axis(1, at = c(3.5, 4, 5, 6, 6.5), line = 2,
                  tick = T, labels = rep("", 5), lwd = 3, lwd.ticks = 0,
                  col = "green")
                barplot(bics, names = c("1", "2"), ylab = "BIC",
                  pch = 16, xlab = "Number of components", log = "y")
                plot(runif(length(chr_exp), 0, 100), chr_exp,
                  pch = 16, col = cellcolor, ylab = "Expression z-score",
                  ylim = c(min(chr_exp), (max(chr_exp) + 2)),
                  xlab = "Cells")
                legend("topright", col = c("black", "red"), c("Non-malignant",
                  "Malignant"), bty = "o", box.col = "darkgreen",
                  cex = 0.65, pch = 16, title = "Pred. via transcript.")
            }
            else {
                if (vis == T) {
                  plot(runif(length(chr_exp), 0, 100), chr_exp,
                    pch = 16, ylab = "Expression z-score", ylim = c(min(chr_exp),
                      (max(chr_exp) + 2)), xlab = "Cells")
                  barplot(bics, names = c("1", "2"), ylab = "BIC",
                    pch = 16, xlab = "Number of components",
                    log = "y")
                  hist(bestmix$posterior[, 1], main = "Posterior probablility distribution\n component 1",
                    xlab = "Posterior probability", breaks = 20,
                    xlim = c(0, 1))
                }
            }
            return(bestmix)
        }
    }
}


plotAll_ensembl = function (mat, normFactor, regions, gene_pos, fname, normal = NULL,
    tumor = NULL, postProb = 0.8, repetitions = 4)
{
    pdf(paste(fname, "_CNVs.pdf", sep = ""))
    loglik = c()
    bic = c()
    lrt = c()
    l = c()
    for (i in 1:nrow(regions)) {
        mixmdl = plotChrEnrichment_ensembl(mat, regions[i, 1], normFactor,
            gene_pos, nrow(regions), normal, tumor, regions[i,
                2], regions[i, 3], postProb = postProb, repetitions = repetitions)
        if (!is.null(mixmdl)) {
            loglik = c(loglik, mixmdl$loglik)
            bic = c(bic, mixmdl$BIC)
            lrt = c(lrt, mixmdl$lrt)
            names(loglik)[length(loglik)] = rownames(regions)[i]
            names(bic)[length(bic)] = paste(rownames(regions)[i],
                "1_comp", sep = "_")
            names(bic)[length(bic) - 1] = paste(rownames(regions)[i],
                "2_comp", sep = "_")
            if (mixmdl$mu[1] > mixmdl$mu[2]) {
                r = mixmdl$posterior[, 1]
            }
            else {
                r = mixmdl$posterior[, 2]
            }
            l = cbind(l, r)
            colnames(l)[ncol(l)] = rownames(regions)[i]
        }
    }
    par(mfrow = c(1, 1))
    barplot(sort(loglik), names = names(sort(loglik)), cex.axis = 0.8,
        cex.names = 0.7, las = 2, ylab = "log-likelihood")
    dev.off()
    bicLRmat = matrix(ncol = 4, nrow = length(loglik))
    bicLRmat[, 1] = bic[seq(1, (length(bic) - 1), 2)]
    bicLRmat[, 2] = bic[seq(2, length(bic), 2)]
    bicLRmat[, 3] = bicLRmat[, 1] - bicLRmat[, 2]
    bicLRmat[, ncol(bicLRmat)] = lrt
    colnames(bicLRmat) = c("BIC 1 component", "BIC 2 components",
        "BIC difference", "LRT adj. p-val")
    rownames(bicLRmat) = names(loglik)
    rownames(l) = colnames(mat)
    write.table(bicLRmat, paste(fname, "BIC_LR.txt", sep = "_"),
        sep = "\t")
    return(l)
}
"""
#suva_expr = as.matrix(read.table("OG_processed_data_portal.txt",sep="\t",header=T,row.names=1,check.names=F))
suva_expr = as.matrix(read.table("MN1_OG_processed_data_portal.txt",sep="\t",header=T,row.names=1,check.names=F))
suva_expr [which(is.na(suva_expr ))]=0
dim(suva_expr)
suva_expr[1:5,1:5]

# Load cluster information
clusters = read.csv('MN1_Cluster_Identities.csv')
head(clusters)
clusters1 = clusters$leiden

#patients=unlist(strsplit(colnames(suva_expr),"_",fixed=TRUE))[seq(1,(3*ncol(suva_expr))-1,3)]
#unique(patients)
#patients[which(patients=="93")]="MGH93"
#patients[which(patients=="97")]="MGH97"

regions=read.table("chromosome_arm_positions_grch38.txt",sep="\t",row.names = 1,header = T)
head(regions,n=5)

# Preprocessing
gene_pos=getGenePositions(rownames(suva_expr), filter = "ensembl_gene_id")
#gene_pos=getGenePositions_ensembl(rownames(suva_expr)
#suva_expr=filterMatrix(suva_expr,gene_pos[,"hgnc_symbol"],minCells=5)
suva_expr=filterMatrix(suva_expr,gene_pos[,"ensembl_gene_id"],minCells=5)
normFactor=calcNormFactors(suva_expr)
l=plotAll(suva_expr,normFactor,regions,gene_pos, gene_pos_index = 1, "MN1_SUVA_CNVs")
# save out, read into scanpy, push into an obs and paint
write.csv(l, "MN1_posterior_probabilities_CONICSmap.csv")
# python load into pandas as dataframe
# obs scanpy, axis=1

hi=plotHistogram(l,suva_expr,clusters=1,zscoreThreshold=4,clusters1)
# 20q specific to cells - paint that onto UMAPs
dev.off()


vg=detectVarGenes(suva_expr,500)

ts=calculateTsne(suva_expr,vg)
plotTsneGene(ts,suva_expr,c("MBP","CSF1R","ALDOC","OLIG1"))

plotTsneProbabilities(ts,suva_expr,l[,"1p"],"Tsne Chr1p")

lrbic=read.table("SUVA_CNVs_BIC_LR.txt",sep="\t",header=T,row.names=1,check.names=F)
colnames(lrbic)
candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>200 & lrbic[,"LRT adj. p-val"]<0.01)]
hi=plotHistogram(l[,candRegions],suva_expr,clusters=4,zscoreThreshold=4,patients)
