# Generate snakemake config file containing:
# fastq_samples:
#    sample_name: sample_path

import os

sample_lst = os.listdir('data/rnaseq/fastq')

sample_lst = [s.replace(".fastq.gz", "") + ": data/rnaseq/fastq/"+s for s in
              sample_lst]

with open("config.yaml", "w") as f:
    f.write("samples:\n  ")
    f.write("\n  ".join(sample_lst))
    f.write("\n")
