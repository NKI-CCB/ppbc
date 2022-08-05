#cd /DATA/share/postpartumbc/data/external/index
#wget http://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#wget http://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz

#gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#gunzip Homo_sapiens.GRCh38.94.gtf.gz

#mkdir STAR_grch38_index


STAR --runThreadN 32 \
--runMode genomeGenerate \
--genomeDir /DATA/share/postpartumbc/data/external/index/STAR_grch38_index \
 --genomeFastaFiles /DATA/share/postpartumbc/data/external/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --sjdbGTFfile /DATA/share/postpartumbc/data/external/index/Homo_sapiens.GRCh38.94.gtf \
 --sjdbOverhang 50 #Read length is 51, ideal value is length - 1
