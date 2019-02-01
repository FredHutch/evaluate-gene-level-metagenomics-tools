workflow makeUniqueFastaHeaders_wf {

  File fasta_input
  
  call makeUniqueFastaHeaders {
    input:
      fasta_input=fasta_input
  }

  output {
    File fasta_output=makeUniqueFastaHeaders.fasta_output
  }

}

task makeUniqueFastaHeaders {

  File fasta_input
  String sample_base = basename(fasta_input, ".gz")

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
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

with open("TEMP", "wt") as fo:
    for header, seq in SimpleFastaParser(gzip.open("${fasta_input}", "rt")):
        header = make_new_name(header)
        fo.write(">" + header + "\n" + seq + "\n")
        written_genes.add(header)

END
    mv TEMP "${sample_base}"
    gzip "${sample_base}"

  }
  output {
    File fasta_output = "${sample_base}.gz"
  }
}
