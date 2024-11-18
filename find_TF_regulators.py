###################################################
## Functions to write PBN networks               ##
##  ______     ______     __  __                 ##
## /\  __ \   /\  ___\   /\ \/\ \                ##
## \ \  __ \  \ \___  \  \ \ \_\ \               ##
##  \ \_\ \_\  \/\_____\  \ \_____\              ##
##   \/_/\/_/   \/_____/   \/_____/              ##
## @Developed by: Plaisier Lab                   ##
##   (https://plaisierlab.engineering.asu.edu/)  ##
##   Arizona State University                    ##
##   242 ISTB1, 550 E Orange St                  ##
##   Tempe, AZ  85281                            ##
## @Author:  Chris Plaisier                      ##
## @License:  GNU GPLv3                          ##
###################################################

from optparse import OptionParser
import sys
import gzip

parser = OptionParser()
usage = "usage: %prog --input FILE --output FILE"
parser = OptionParser(usage=usage)
parser.add_option('-i', '--input', dest='input',
                  help='input gene list (symbols)', metavar='FILE')
parser.add_option('-o', '--output', dest='output',
                  help='file to dump output', metavar='FILE')

(options, args) = parser.parse_args()

if not options.input or not options.output:
    parser.print_help()
    sys.exit(1)

from scipy.stats import hypergeom
import json

def enrichment(geneList1, geneList2, allGenes):
    tmp = set(geneList1).intersection(geneList2)
    x = len(tmp)
    M = allGenes
    n = len(geneList1)
    N = len(geneList2)
    return [x,tmp,M,n,N,1-hypergeom.cdf(x, M, n, N)]

# Read in gene symbol conversion to Entrez IDs
symbol2entrezId = {}
entrezId2symbol ={}
with open('libs/mart_export.txt','r') as inFile:
    inFile.readline() # Get rid of header
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split('\t')
        if len(splitUp)==3 and not (splitUp[1]=='' or splitUp[2]==''):
            symbol2entrezId[splitUp[1]] = splitUp[2]
            entrezId2symbol[splitUp[2]] = splitUp[1]

# Read in humanTFs_All.csv with <TF Name>,<Entrez ID>
tfName2entrezId = {}
with open('libs/humanTFs_All.csv','r') as inFile:
    inFile.readline() # Get rid of header
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        tfName2entrezId[splitUp[0]] = splitUp[2]

# Read in tfFamilies.csv for expanded list of possible TFs
tfFamilies = {}
with open('libs/tfFamilies.csv','r') as inFile:
    inFile.readline() # Get rid of header
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(',')
        tmp = splitUp[2].split(' ')
        tfFamilies[splitUp[0]] = tmp

# Load up the TF target genes
allGenes = set()
tfTargets = {}
with gzip.open('libs/tfbsDb_plus_and_minus_5000_entrez.json.gz','rb') as inFile:
    tfTargets = json.load(inFile)
allGenes = set([item for sublist in tfTargets.values() for item in sublist])

# Load up gene list for comparison
geneList1 = []
with open(options.input,'r') as inFile:
    geneList1 = [symbol2entrezId[i.strip()] for i in inFile.readlines() if i.strip() in symbol2entrezId]

# Do enrichment analysis
output = ['TF,TF Symbol,TF Family,TF Family Members,Overlap,Overlapping,All Genes,Mapped Input Gene List,TF Targets,P-Value']
for tf1 in tfTargets:
    res1 = enrichment(geneList1, tfTargets[tf1], len(allGenes))
    tfName = ''
    tfFam = ''
    posNames = ''
    if tf1 in tfName2entrezId:
        tfName = entrezId2symbol[tfName2entrezId[tf1]]
        for fam in tfFamilies:
            if tfName2entrezId[tf1] in tfFamilies[fam]:
                tfFam = fam
        if not tfFam=='':
            posNames = ' '.join([entrezId2symbol[i] for i in tfFamilies[tfFam] if i in entrezId2symbol])
    output.append(','.join([tf1, tfName, tfFam, posNames, str(res1[0]), ' '.join([entrezId2symbol[i] for i in res1[1] if i in entrezId2symbol]), str(res1[2]), str(res1[3]), str(res1[4]), str(res1[5])]))
with open(options.output,'w') as outFile:
    outFile.write('\n'.join(output))

