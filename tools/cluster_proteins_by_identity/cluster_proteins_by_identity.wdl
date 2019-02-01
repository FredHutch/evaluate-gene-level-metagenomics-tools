workflow clusterProteinsAbunds {

  File ref_fasta_in
  File ref_abund_in
  File detected_fasta_in
  String identity="0.9"
  String overlap="0.9"

  call clusterProteins as cluster_ref {
    input:
      fasta_in=ref_fasta_in,
      identity=identity,
      overlap=overlap
  }

  call clusterProteins as cluster_det {
    input: 
      fasta_in=detected_fasta_in,
      identity=identity,
      overlap=overlap
  }

  call clusterAbund {
    input:
      abund_in=ref_abund_in,
      groups=cluster_ref.tsv_out,
      identity=identity
  }

  output {
    File ref_fasta_out=cluster_ref.fasta_out
    File ref_abund_out=clusterAbund.abund_out
    File detected_fasta_out=cluster_det.fasta_out
  }

}

task clusterProteins {
  File fasta_in
  String identity="0.9"
  String overlap="0.9"
  String base=basename(fasta_in)
  String cpu="8"

  runtime {
    docker: "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    memory: "4G"
    cpu: cpu
  }

  command {
    set -e;

    mmseqs createdb ${fasta_in} DB;
    mmseqs cluster DB clu tmp --min-seq-id ${identity} --max-seqs 1000000 -c ${overlap} --threads ${cpu};
    mmseqs createtsv DB DB clu ${base}.clusters.tsv;
    mmseqs result2repseq DB clu ${base}.rep;
    mmseqs result2flat DB DB ${base}.rep ${base}.rep.fasta --use-fasta-header

  }
  output {
    File fasta_out = "${base}.rep.fasta"
    File tsv_out = "${base}.clusters.tsv"
  }

}

task clusterAbund {
  File abund_in
  File groups
  String identity
  String base=basename(abund_in, ".csv")

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    python << END

from collections import defaultdict

# Read in the cluster membership
clusters = dict([
  line.rstrip("\n").split("\t")[:2]
  for line in open("${groups}", "rt")
])

# Add up the abundances for each cluster
clust_abund = defaultdict(float)
for line in open("${abund_in}", "rt"):
    if ',' not in line:
        continue
    member, abund = line.rstrip("\n").split(",", 2)[:2]
    assert member in clusters
    clust_abund[clusters[member]] += float(abund)

with open("${base}.clust.${identity}.csv", "wt") as fo:
    for k, v in clust_abund.items():
        fo.write(k + ',' + str(v) + '\n')

END
  }

  output {
    File abund_out = "${base}.clust.${identity}.csv"
  }

}