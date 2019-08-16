# Generate symlinks to fastq files
import os

f_paths = open('/DATA/share/postpartumbc/data/metadata/sample-paths.txt', 'r')
sample_paths = f_paths.read().splitlines()


sample_list = []
for s in sample_paths:
    src = '/DATA/share/postpartbc/' + s
    dst = "/DATA/share/postpartumbc/data/RAW/fastq_dec2018/" + s.split('/')[-1]
    os.symlink(src, dst)
