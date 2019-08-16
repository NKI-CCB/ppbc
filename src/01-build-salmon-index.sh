# Build Salmon index from Ensembl release 94
# Last updated: Tue Nov 26 2018, Gergana Bounova
# Adapted from Gergana Bounova
# Index rebuilt under version 0.14 of salmon in jun2019 upon receipt of additional samples

cd ~/projects/postpartbc/data/external/index

# Making an transcriptome index:
# Downloaded ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# and ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
# and (unzipped) concatenated (and zipped) them to
# Homo_sapiens.GRCh38.cdna.all.ncrna.fa.

wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.ncrna.fa.gz

cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > Homo_sapiens.GRCh38.cdna.all.ncrna.fa

rm Homo_sapiens.GRCh38.cdna.all.fa
rm Homo_sapiens.GRCh38.ncrna.fa

# Build the index
salmon index -t Homo_sapiens.GRCh38.cdna.all.ncrna.fa -i grch38_index --type quasi -k 23
