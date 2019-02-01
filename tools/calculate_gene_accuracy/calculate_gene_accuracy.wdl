import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/1288836ffeb21f3939249e712b59774567ffc11c/WDL/align-proteins-diamond/align-proteins-diamond.wdl" as dmd
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/f86bca3613d96bae828318d9bedbbe4e4716670e/tools/make_unique_fasta_headers/make_unique_fasta_headers.wdl?token=AE-VSOIyw6rQ7HliYuS5j5AODY17nmxgks5cXeWfwA" as unique_fasta

workflow calculateGeneAccuracy {

  File ref_fasta
  File ref_abund
  File detected_fasta
  String method_label="dummy_method"

  call dmd.MakeDiamondDatabase {
    input:
      fasta=ref_fasta
  }

  call unique_fasta.makeUniqueFastaHeaders {
    input:
      fasta_input=detected_fasta
  }

  call dmd.RunDiamond {
    input:
      db=MakeDiamondDatabase.db,
      query=makeUniqueFastaHeaders.fasta_output
  }

  call calculateGeneAccuracyFromAln {
    input:
      detected_fasta=makeUniqueFastaHeaders.fasta_output,
      aln=RunDiamond.aln,
      ref_abund=ref_abund,
      method_label=method_label
  }
  
  output {
    File accuracy=calculateGeneAccuracyFromAln.accuracy
  }

}

task calculateGeneAccuracyFromAln {
  File detected_fasta
  File aln
  File ref_abund
  String method_label
  String output_name=basename(detected_fasta, ".fasta.gz")

  runtime {
    docker: "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    memory: "2G"
    cpu: "1"
  }

  command {
    set -e; 
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

  }
  output {
    File accuracy = "${output_name}.accuracy.tsv"
  }

}