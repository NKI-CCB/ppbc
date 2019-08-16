# Generate symlinks to fastq files
import os

f_paths = open('/DATA/share/postpartumbc/data/metadata/sample_paths_jun2019.txt', 'r')
sample_paths = f_paths.read().splitlines()


sample_list = []
for s in sample_paths:
    src = '/DATA/share/postpartbc/fastq_jun2019/' + s
    dst = "/DATA/share/postpartumbc/data/RAW/fastq_jun2019/" + s.split('/')[-1]
    os.symlink(src, dst)
