#!/bin/bash

set -e

python << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
import gzip

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

read_counts = defaultdict(int)
read_assignments = dict()

for line in gzip_safe_open("${aln}"):
    line = line.rstrip("\n").split("\t")
    if len(line) < 5:
        continue
    read_assignments[line[0]] = line[1]
    read_counts[line[0]] += 1

detected_genes = set([
  v
  for k, v in read_assignments.items()
  if read_counts[k] == 1
])

written_genes = dict()
for header, seq in SimpleFastaParser(gzip_safe_open("${fasta_in}")):
    if header in detected_genes:
        written_genes[header] = seq

assert len(written_genes) == len(detected_genes)

with open("${sample_base}.fasta", "wt") as fo:
    for k, v in written_genes.items():
        fo.write(">" + k + "\n" + v + "\n")

END

gzip "${sample_base}.fasta"
