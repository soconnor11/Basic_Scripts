# file dir = plaisierlab1/home/cplaisier/ccAF_newData/usftp21.novogene.com/raw_data

import os.path
from subprocess import *

# Sample information
samples = [ 'hBMMSC_CKDL210009543-1a-SI_GA_G1_H57KLDSX2', 'hESC_CKDL210009544-1a-SI_GA_H2_H57KLDSX2' ]

# Run cellranger count pipeline for each sample
for sample1 in samples:
    if not os.path.exists('/files/cellranger_output/'+sample1):
        cmd = 'cellranger count --id='+sample1+' --transcriptome=/opt/refdata-gex-GRCh38-2020-A --fastqs=/files --sample='+sample1+' --expect-cells=4000 --localcores=16'
        print('Running = '+cmd)
        with open('stdout.txt','w') as stdoutFile, open('stderr.txt','w') as stderrFile:
            cellranger = Popen(cmd, shell=True, stdout=stdoutFile, stderr=stderrFile)
            cellranger.communicate()
