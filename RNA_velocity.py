
#------------------------------------------------------
# RNA velocity - 08/05/21
#---------------------------------------------------

#------------------------------------------------------
# Set up post-processing bam files for velocyto
#---------------------------------------------------

docker run -it -v '/home/soconnor/ccAF_newData:/files/hBMMSC_CKDL210009543-1a-SI_GA_G1_H57KLDSX2/outs/' cplaisier/scrna_seq_velocity_6_7_2021
samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam  # takes a while


#------------------------------------------------------
# Run velocyto - dump out loom file
#---------------------------------------------------

docker run -it -v '/home/soconnor/ccAF_newData:/files/' cplaisier/scrna_seq_velocity_6_7_2021
velocyto run10x -m GRCh38_rmsk.gtf hBMMSC_CKDL210009543-1a-SI_GA_G1_H57KLDSX2/outs genes_GRCh38.gtf

#ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 88 from C header, got 80 from PyObject



~/anaconda3/bin/velocyto run10x -m ~/Downloads/hg19_rmsk.gtf ~/Dropbox_ASU/Feldman_scRNA_seq/u5_all/U5_all ~/Downloads/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
