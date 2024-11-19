
#------------------------------------------------------
# RNA velocity - 08/05/21
#---------------------------------------------------


#------------------------------------------------------
# Generate loom file from velocyto
#---------------------------------------------------
docker run -it -v '/home/soconnor/ccAF_newData:/files/hBMMSC_CKDL210009543-1a-SI_GA_G1_H57KLDSX2/outs/' cplaisier/scrna_seq_velocity_6_7_2021
samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam  # takes a while

~/anaconda3/bin/velocyto run10x -m ~/ccAF_newData/GRCh38_rmsk.gtf ~/ccAF/hBMMSC_CKDL210009543-1a-SI_GA_G1_H57KLDSX2/outs ~/ccAF_newData/genes_GRCh38.gtf
