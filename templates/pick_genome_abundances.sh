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

# Figure out which column has the organism name
assert "Organism" in header
organism_ix = header.index("Organism")

# Make a list of the fields in each line
list_of_lines = [
    l.rstrip("\\n").split("${sep}")
    for l in f
]

# Get the organism name and RefSeq FTP path for each genome
all_genomes = [
    (l[organism_ix], l[refseq_ftp_ix])
    for l in list_of_lines
    if len(l[refseq_ftp_ix]) > 1 and l[refseq_ftp_ix].startswith('ftp://')
]

print("There are " + str(len(all_genomes)) + " genomes available")
assert num_genomes <= len(all_genomes)

# Randomly shuffle the genomes
random.shuffle(all_genomes)

# Only take a single genome per species
filtered_genomes = []
species_seen = set([])
for org_name, ftp_path in all_genomes:
    species_name = " ".join(org_name.split(" ")[:2]) if len(org_name.split(" ")) > 2 else org_name
    if species_name not in species_seen:
        filtered_genomes.append(ftp_path)
        species_seen.add(species_name)

print("There are " + str(len(filtered_genomes)) + " genomes after removing species duplicates")

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
      filtered_genomes[:num_genomes],
      abund_list
    ):
        print(k, v)
        fo.write(k + ',' + str(v) + '\\n')

END
