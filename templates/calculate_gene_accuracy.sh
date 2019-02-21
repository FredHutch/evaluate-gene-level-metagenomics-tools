#!/bin/bash

set -e

python << END

import gzip
import pandas as pd

# Read in the table with the depth for each reference protein
ref_abund = pd.read_table("${ref_abund}", sep=",", header=None)

# Read in the list of detected genes
detected_genes = set([
    header[1:].split(" ")[0].rstrip("\n")
    for header in gzip.open("${detected_fasta}", "rt")
    if header[0] == ">"
])

# Read in the set of alignments
aln = pd.read_table("${aln}", sep="\t", header=None)

# Make sure there are >0 rows
assert aln.shape[0] > 0

# Make sure there are 12 columns
assert aln.shape[1] == 12

# Sort by alignment score
aln.sort_values(by=11, inplace=True, ascending=False)

# Take the top hit for every query
aln = aln.groupby(0).head(1)

# Format the output as a list of dicts
output = []

# Add the number of false positives (with no alignment)

output.append(dict([
    ("method", "${method_label}"),
    ("depth", "all"),
    ("false_positive", len(detected_genes) - aln[0].isin(detected_genes).sum())
]))

# Add the TP, FN, and DUP at each level of depth
for d, df in ref_abund.groupby(1):
    d = float(d)
    tp = len(set(df[0]) & set(aln[1].tolist()))
    fn = df.shape[0] - tp
    dup = aln[1].isin(df[0]).sum() - tp
    output.append(dict([
        ("method", "${method_label}"),
        ("depth", d),
        ("true_positive", tp),
        ("false_negative", fn),
        ("duplicate", dup)
    ]))

output = pd.DataFrame(output)
output.to_csv("${output_name}.accuracy.tsv", sep="\t", index=None)

END