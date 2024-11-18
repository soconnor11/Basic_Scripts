#docker run -it -v '/home/soconnor/ccAF_newData:/files' cplaisier/sygnal

#------------------------
# Set up / imports
#-----------------------

library(stringr)
library(GOSemSim)
library(getopt)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
hsGO = godata('org.Hs.eg.db',ont='BP')

# read command line arguments
spec = matrix(c('resdir', 'd', 1, 'character',
                'savedir', 's', 1, 'character',
                'savedirFilt', 'f', 1, 'character',
                'help', 'h', 0, 'logical'), byrow=TRUE, ncol=4)
opt <- getopt(spec)
if (is.null(opt$resdir) || is.null(opt$savedir) || is.null(opt$savedirFilt) || !is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

#resdir = 'redo_analysis/TRN/markerEnrichment/filtered/'

#------------------------
# Hallmark enrichment analysis
#-----------------------

# Load reference hallmarks
#l1 = list()
#l1$glial_cell_differentiation = c('GO:0010001')
#l1$cell_quiescence = c('GO:0044838') # NAs for all
#l1$G0_G1_transition = c('GO:0045023')
#l1$G1_S_transition = c('GO:0000082')
#l1$dna_replication = c('GO:0006260')
#l1$G2M_transition = c('GO:0000086')
#l1$cell_division = c('GO:0051301')
#l1$nuclear_envelope_reassembly = c('GO:0031468')

l1 = list()
l1$SelfSufficiencyInGrowthSignals = c('GO:0009967','GO:0030307','GO:0008284','GO:0045787','GO:0007165')
l1$InsensitivityToAntigrowthSignals = c('GO:0009968','GO:0030308','GO:0008285','GO:0045786','GO:0007165')
l1$EvadingApoptosis = c('GO:0043069','GO:0043066')
l1$LimitlessReplicativePotential = c('GO:0001302','GO:0032206','GO:0090398')
l1$SustainedAngiogenesis = c('GO:0045765','GO:0045766','GO:0030949','GO:0001570')
l1$TissueInvasionAndMetastasis = c('GO:0042060','GO:0007162','GO:0033631','GO:0044331','GO:0001837','GO:0016477','GO:0048870','GO:0007155')
l1$GenomeInstabilityAndMutation = c('GO:0051276','GO:0045005','GO:0006281')
l1$TumorPromotingInflammation = c('GO:0002419','GO:0002420','GO:0002857','GO:0002842','GO:0002367','GO:0050776')
l1$ReprogrammingEnergyMetabolism = c('GO:0006096','GO:0071456')
l1$EvadingImmuneDetection = c('GO:0002837','GO:0002418','GO:0002367','GO:0050776')

# Load data
#d1 = read.csv('ccAF_newData/redo_analysis/TRN/markerEnrichment/filtered/G0_tf_analysis_output_for_Hallmarks.csv',header=T,row.names=1)
#file1 = "G0_tf_analysis_output_for_Hallmarks.csv"

files <- list.files(path = opt$resdir, pattern = ".csv")
for(file1 in files) {
  d1 = read.csv(file = paste0(opt$resdir,file1),header=T,row.names=1,sep=',')
  tmp <- file1
  tag <- str_split(tmp, "_")[[1]][1]
  l2 = list()
  for(cluster in rownames(d1)) {
      l2[[cluster]] = intersect(trimws(strsplit(as.character(d1[cluster,1]),',')[[1]]),hsGO@geneAnno$GO)
  }
  hallmarks = matrix(ncol=length(names(l1)), nrow=length(names(l2)), dimnames=list(names(l2), names(l1)))
  for(cluster in names(l2)) {
      if (!(length(l2[[cluster]])==0)) {
          for(hallmark in names(l1)) {
              #print(cluster)
              #d2 = getTermSim(c(l2[[cluster]],l1[[hallmark]]),method='JiangConrath')
              #hallmarks[cluster,hallmark] = max(d2[1:length(l2[[cluster]]),-(1:length(l2[[cluster]]))],na.rm=T)
              hallmarks[cluster,hallmark] = mgoSim(l2[[cluster]], l1[[hallmark]], semData=hsGO, measure='Jiang', combine='max')
          }
      }
  }
  write.csv(hallmarks, file.path(opt$savedir, paste0(tag,"_jiangConrath_hallmarks_cancer.csv")))
  tmp2 = subset(data.frame(hallmarks), rowSums(data.frame(hallmarks)[1:10] >= 0.8) > 0)
  if (!(dim(tmp2)[1[1]]==0)) {
      write.csv(tmp2, file.path(opt$savedirFilt, paste0(tag,"_jiangConrath_hallmarks_cancer.csv")))
      }
}
