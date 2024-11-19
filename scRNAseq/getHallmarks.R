#docker run -it -v '/home/soconnor/ccAF_newData:/files' cplaisier/sygnal

#------------------------
# Set up / imports
#-----------------------

library(stringr)
suppressMessages(library(GOSemSim))
suppressMessages(library(getopt))
​hsGO = godata('org.Hs.eg.db',ont='BP')

# read command line arguments
spec = matrix(c(
  ‘resdir’, ‘d’, 1, ‘character’,
  ‘help’, ‘h’, 0, ‘logical’
), byrow=TRUE, ncol=4)
opt <- getopt(spec)
if (is.null(opt$resdir) || !is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

#resdir = 'ccAF_newData/redo_analysis/TRN/markerEnrichment/filtered/'

#------------------------
# Hallmark enrichment analysis
#-----------------------

# Load reference hallmarks
l1 = list()
l1$glial_cell_differentiation = c('GO:0010001')
l1$G0_G1_transition = c('GO:0045023')
l1$G1_S_transition = c('GO:0000082')
l1$dna_replication = c('GO:0006260')
l1$G2M_transition = c('GO:0000086')
l1$cell_division = c('GO:0051301')
l1$nuclear_envelope_reassembly = c('GO:0031468')
​

# Load data
#d1 = read.csv('redo_analysis/markerEnrichment/S_hallmark_analysis_output_filtered.csv',header=T,row.names=1)

files <- list.files(path = opt$resdir, pattern = ".csv")
for(file1 in files) {
  d1 = read.csv(file = paste0(opt$resdir,file1),header=T,row.names=1,sep=',')
  tmp <- file1
  tag <- str_split(tmp, "_")[[1]][1]
  l2 = list()
  for(cluster in rownames(d1)) {
      l2[[cluster]] = intersect(strsplit(as.character(d1[cluster,1]),',')[[1]],hsGO@geneAnno$GO)
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
  write.csv(hallmarks, file.path(opt$resdir, paste(tag,"_jiangConrath_hallmarks.csv")))
}
