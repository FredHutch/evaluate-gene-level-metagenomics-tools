#!/bin/bash

set -e;

tar xvf "${genome_tar}"

python << END
import os
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Write out the gene sequences and the simulated depths
fastq = gzip.open("genes.${random_seed}.fastp.gz", "wt")
abund = gzip.open("genes_abund.${random_seed}.csv.gz", "wt")

# Read in the depth of sequencing for each genome
for line in open("${genome_abund_csv}", "rt"):
    if "," not in line:
        continue

    url, depth = line.rstrip("\\n").split(",", 1)

    depth = float(depth)
    acc = url.split("/")[-1]

    genes_fp = acc + ".faa.gz"
    assert os.path.exists(genes_fp)

    for header, seq in SimpleFastaParser(gzip.open(genes_fp, "rt")):
        header = header.split(" ")[0].split("\t")[0]
        fastq.write(">" + header + '\\n' + seq + '\\n')
        abund.write(str(header) + ',' + str(depth) + '\\n')

fastq.close()
abund.close()

END
