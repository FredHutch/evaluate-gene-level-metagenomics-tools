workflow extractDetectedGenesFAMLI_wf {

  File json_input
  File fasta_input
  
  call extractDetectedGenesFAMLI {
    input:
      json_input=json_input,
      fasta_input=fasta_input
  }

  output {
    File fasta_output=extractDetectedGenesFAMLI.fasta_output
  }

}

task extractDetectedGenesFAMLI {

  File json_input
  File fasta_input
  String sample_base = basename(json_input, ".json.gz")

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    python3 << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import json

# Read in all of the sequence records from the FAMLI output
detected_genes = set([
    r["id"]
    for r in json.load(gzip.open("${json_input}", "rt"))
])

written_genes = set([])

with open("${sample_base}.fasta", "wt") as fo:
    for header, seq in SimpleFastaParser(gzip.open("${fasta_input}", "rt")):
        header = header.split(" ")[0]
        if header in detected_genes:
            assert header not in written_genes, "Duplicate gene names"
            fo.write(">" + header + "\n" + seq + "\n")
            written_genes.add(header)

assert len(written_genes) == len(detected_genes), "Not all genes were found in the input"

END

    gzip ${sample_base}.fasta

  }
  output {
    File fasta_output = "${sample_base}.fasta.gz"
  }
}
