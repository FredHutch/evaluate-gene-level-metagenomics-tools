#!/bin/bash

set -e

python3 << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip

written_genes = set([])

def make_new_name(n):
    n = n.split(" ", 1)[0].split("\t", 1)[0]
    if n not in written_genes:
        return n
    ix = 0
    while n + "_" + str(ix) in written_genes:
        ix += 1
    return n + "_" + str(ix)

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

ix = 0
with open("${fasta_input}.clean.fasta", "wt") as fo:
    for header, seq in SimpleFastaParser(gzip_safe_open("${fasta_input}", "rt")):
        ix += 1
        header = make_new_name(header)
        fo.write(">" + header + "\\n" + seq + "\\n")
        written_genes.add(header)
assert ix > 0, "No sequences were found"

END

gzip "${fasta_input}.clean.fasta"
