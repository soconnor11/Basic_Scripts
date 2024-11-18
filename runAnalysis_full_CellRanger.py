# file dir = plaisierlab1/home/cplaisier/ccAF_newData/usftp21.novogene.com/raw_data
# docker: cplaisier/cellranger_7.0.1


#docker run -it -v '/media/omics1/ccNN/usftp21.novogene.com/01.RawData/:/files' cplaisier/cellranger_7.0.1
docker run -it -v '/media/Qiu/SkMSC/usftp21.novogene.com/01.RawData/:/files' cplaisier/cellranger_7.0.1


import os.path
from subprocess import *

# check md5sums under data directory
#md5sum -c MD5.txt

# Sample information
#samples = { 'ATCC_hBMMSC': 'ATCC_hBMMSC_CKDL220021127-1A_HF7KHDSX5', 'MCF10A': 'MCF10A_CKDL220021125-1A_HF7KHDSX5', 'ScienCell_hBMMSC': 'ScienCell_hBMMSC_CKDL220021126-1A_HFC23DSX5' }
samples = { 'HSkMSC_cr': 'HSkMSC_CKDL240007875-1A_H2T5KDSXC' }
# Make sure the folder holding all of the fastq have the same name

# Run cellranger count pipeline for each sample
for sample1 in samples:
    if not os.path.exists('/files/cellranger_output/'+sample1):
        cmd = 'cellranger count --id='+sample1+' --transcriptome=/opt/refdata-gex-GRCh38-2020-A --fastqs=/files --sample='+samples[sample1]+' --expect-cells=4000 --localcores=16'
        #cmd = 'cellranger count --id='+sample1+' --transcriptome=/opt/refdata-gex-GRCh38-2020-A --fastqs=/opt/HSkMSC --sample='+samples[sample1]+' --expect-cells=4000 --localcores=16'
        print('Running = '+cmd)
        with open('stdout.txt','w') as stdoutFile, open('stderr.txt','w') as stderrFile:
            cellranger = Popen(cmd, shell=True, stdout=stdoutFile, stderr=stderrFile)
            cellranger.communicate()
