#!/bin/bash

set -e

python3 << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import json

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

# Read in all of the sequence records from the FAMLI output
detected_genes = set([
    r["id"]
    for r in json.load(gzip_safe_open("${json_input}"))
])

written_genes = set([])

with open("${json_input}.faa", "wt") as fo:
    for header, seq in SimpleFastaParser(gzip_safe_open("${fasta_input}")):
        header = header.split(" ")[0]
        if header in detected_genes:
            assert header not in written_genes, "Duplicate gene names"
            fo.write(">" + header + "\\n" + seq + "\\n")
            written_genes.add(header)

assert len(written_genes) == len(detected_genes), "Not all genes were found in the input"

END

(( \$(cat "${json_input}.faa" | wc -l) > 1 ))

gzip "${json_input}.faa"
