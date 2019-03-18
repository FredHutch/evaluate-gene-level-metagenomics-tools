#!/bin/bash

set -e

python3 << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Seq
import gzip
import json
import os

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

# Get the detected genes
detected_genes = set([
    line.split("\\t")[0].split("|")[0]
    for line in gzip_safe_open("${gene_tsv}", "rt")
    if "\\t" in line and "UniRef" in line and "unknown" not in line
])
print(" ".join(["Detected", str(len(detected_genes)), "genes from the HUMAnN2 output"]))
assert len(detected_genes) > 0, "No genes detected in the input"

# Write out the detected genes
written_genes = set([])
with open("humann2.detected.fasta", "wt") as fo:

    # Look through the UniRef database
    with gzip_safe_open("${ref_fasta}", "rt") as fi:
        for header, seq in SimpleFastaParser(fi):
            header = header.split(" ", 1)[0].split("|", 1)[0]
            if header in detected_genes:
                written_genes.add(header)
                fo.write(">" + header + "\\n" + seq + "\\n")

    print(" ".join(["Wrote out", str(len(written_genes)), "genes in the UniRef database"]))

    print("Looking through the chocophlan database")

    # Look through the Chocophlan database
    for fname in os.listdir("chocophlan"):

        # Skip files that aren't genome references
        if fname.endswith(".ffn.gz") is False:
            continue

        # Skip files once we've retrieved all of the detected genes
        if len(written_genes) == len(detected_genes):
            continue

        # Iterate over every item in each reference file
        for header, seq in SimpleFastaParser(gzip_safe_open("chocophlan/" + fname)):

            # Split up the header by '|' to find the reference name used in the output
            header = [f for f in header.split("|") if f.startswith("UniRef")]
            if len(header) == 0:
                continue
            # Take the first field
            header = header[0]

            # Check if we need this particular reference
            if header in detected_genes and header not in written_genes:
                written_genes.add(header)
                fo.write(">" + header + "\\n" + Seq.translate(seq, table=11) + "\\n")

print(" ".join(["Wrote out", str(len(written_genes)), "genes overall"]))

assert len(written_genes) == len(detected_genes), (len(written_genes), len(detected_genes))

END
