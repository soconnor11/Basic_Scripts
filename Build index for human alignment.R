
# Build index for human alignment
mkdir /GRCh38.p13
cd /GRCh38.p13
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
pigz -d *.gz
mkdir /index
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /index --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v32.primary_assembly.annotation.gtf --sjdbOverhang 100

cd /

# Get rid of sequence information to reduce image size
rm -rf /GRCh38.p13/GRCh38.primary_assembly.genome.fa
