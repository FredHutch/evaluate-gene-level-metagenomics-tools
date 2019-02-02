workflow checkUniqueFastaHeaders_wf {

  File fasta_input
  
  call checkUniqueFastaHeaders {
    input:
      fasta_input=fasta_input
  }

  output {
    File fasta_output=checkUniqueFastaHeaders.fasta_output
  }

}

task checkUniqueFastaHeaders {

  File fasta_input
  String base = basename(fasta_input)

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    python3 << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip

gene_names = set([])

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

ix = 0
for header, seq in SimpleFastaParser(gzip_safe_open("${fasta_input}", "rt")):
    ix += 1
    assert header not in gene_names
    gene_names.add(header)
assert ix > 0, "No sequences were found"

END

  cp "${fasta_input}" "${base}"

  }
  output {
    File fasta_output = "${base}"
  }
}
