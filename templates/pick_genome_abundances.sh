#!/bin/bash


set -e;

echo "Simulating index ${ix}"

python << END
import gzip
import random
import numpy as np

num_genomes = int(${num_genomes})
mean_depth = float(${mean_depth})
assert mean_depth > 0.01 and mean_depth < 100, "Please pick a reasonable depth"

# Open a file handle
f = gzip.open("${genome_list}", "rt")

# Read in the header
header = f.readline().rstrip("\\n").split("${sep}")

# Figure out which column has the RefSeq accession
assert "RefSeq FTP" in header
refseq_ftp_ix = header.index("RefSeq FTP")

# Get the complete list of RefSeq FTP paths
all_genomes = [
  l.rstrip("\\n").split("${sep}")[refseq_ftp_ix]
  for l in f
]
all_genomes = [g for g in all_genomes if len(g) > 1 and g.startswith('ftp://')]

print("There are " + str(len(all_genomes)) + " genomes available")
assert num_genomes <= len(all_genomes)

# Randomly shuffle the genomes
random.shuffle(all_genomes)

# Randomly pick a set of abundances
log_mean = np.log10(mean_depth)
abund_list = np.random.normal(
  loc=log_mean,
  scale=${log_std},
  size=num_genomes
)

# Transform back to linear space
abund_list = [
  min(10**x, ${max_depth})
  for x in abund_list
]

# Write out to a file
with open("genome_abund.${ix}.csv", "wt") as fo:
    for k, v in zip(
      all_genomes[:num_genomes],
      abund_list
    ):
        print(k, v)
        fo.write(k + ',' + str(v) + '\\n')

END
