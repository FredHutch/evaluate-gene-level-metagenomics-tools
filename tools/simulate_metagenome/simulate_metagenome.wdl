import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/be533c85a18a7807313f84518f80a857bdfd6fd4/tools/make_unique_fasta_headers/make_unique_fasta_headers.wdl?token=AE-VSOKDLxrQh5a8n42Do3DhQNFY5tREks5cXh9LwA" as clean_headers

workflow SimulateMetagenome {

  File genome_list
  Int num_metagenomes = 5
  String num_genomes = "100"
  String mean_depth = "1"
  String max_depth = "20"
  
  scatter (ix in range(num_metagenomes)){

    call PickGenomeAbundances {
      input:
        genome_list=genome_list,
        num_genomes=num_genomes,
        mean_depth=mean_depth,
        max_depth=max_depth,
        random_seed=ix
    }

    call DownloadGenomes {
      input:
        genome_abund_csv=PickGenomeAbundances.genome_abund_csv
    }

    call SimulateGenomes {
      input:
        genome_abund_csv=PickGenomeAbundances.genome_abund_csv,
        genome_tar=DownloadGenomes.genome_tar
    }

    call InterleaveFastqs {
      input:
        all_reads_tar=SimulateGenomes.all_reads_tar
    }

    call MakeGeneAbundances {
      input:
        genome_abund_csv=PickGenomeAbundances.genome_abund_csv,
        genome_tar=DownloadGenomes.genome_tar
    }

  }

  output {
    Array[File] abund_csv=MakeGeneAbundances.gene_abund_csv
    Array[File] genes_fastp=MakeGeneAbundances.genes_fastp
    Array[File] reads_fastq=InterleaveFastqs.reads_fastq
  }

}

task MakeGeneAbundances {

  File genome_abund_csv
  File genome_tar
  
  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    set -e;

    tar xvf "${genome_tar}"

    python << END
import os
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Write out the gene sequences and the simulated depths
fastq = gzip.open("genes.fastp.gz", "wt")
abund = gzip.open("genes_abund.csv.gz", "wt")

# Read in the depth of sequencing for each genome
for line in open("${genome_abund_csv}", "rt"):
    if "," not in line:
        continue

    url, depth = line.rstrip("\n").split(",", 1)

    depth = float(depth)
    acc = url.split("/")[-1]

    genes_fp = acc + ".faa.gz"
    assert os.path.exists(genes_fp)

    for header, seq in SimpleFastaParser(gzip.open(genes_fp, "rt")):
        header = header.split(" ")[0].split("\t")[0]
        fastq.write(">" + header + '\n' + seq + '\n')
        abund.write(str(header) + ',' + str(depth) + '\n')

fastq.close()
abund.close()

END

  }

  output {
    File gene_abund_csv="genes_abund.csv.gz"
    File genes_fastp="genes.fastp.gz"
  }

}

task SimulateGenomes {
  
  File genome_abund_csv
  File genome_tar
  Int read_length = 100
  Int mflen = 300
  Int sdev = 75
  
  runtime {
    docker: "quay.io/biocontainers/art@sha256:14f44c1cf099f6b55922aaa177c926c993733dba75ef4eb4dcf53442e3b5f96e"
    memory: "4G"
    cpu: "1"
  }

  command {
    set -e;

    tar xvf "${genome_tar}"

    cat "${genome_abund_csv}" | tr ',' '\t' | while read url depth; do
        acc="$(echo $url | sed 's/.*\///')"
        [[ -s $acc.fna.gz ]]

        gunzip $acc.fna.gz

        art_illumina --fcov $depth --in $acc.fna --len ${read_length} --mflen ${mflen} --sdev ${sdev} --noALN --out $acc. --paired --quiet
        gzip *[12].fq

    done

    tar cvf all_reads.tar *.fq.gz

  }

  output {
    File all_reads_tar="all_reads.tar"
  }

}

task InterleaveFastqs {

  File all_reads_tar

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    set -e;

    tar xvf "${all_reads_tar}"

    for r1 in *1.fq.gz; do

        r2="$(echo $r1 | sed 's/1.fq.gz/2.fq.gz/')"

        if [[ -s $r2 ]]; then

            python << END
import gzip

with gzip.open("$r1", "rt") as f1, gzip.open("$r2", "rt") as f2, open("reads.fastq", "wt") as fo:
    while True:
        line = f1.readline()
        if line.strip() == "":
            break
        fo.write(line)
        
        for i in range(3):
            fo.write(f1.readline())
        
        for i in range(4):
            fo.write(f2.readline())
END

        rm "$r1" "$r2"

        fi

    done

  gzip reads.fastq

  }

  output {
    File reads_fastq = "reads.fastq.gz"
  }

}

task DownloadGenomes {

  File genome_abund_csv

  runtime {
    docker: "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    memory: "4G"
    cpu: "1"
  }

  command {
    set -e;

    ending="_genomic.fna.gz"
    unzip_ending="_genomic.fna"

    cat ${genome_abund_csv} | sed 's/,.*//' | while read url; do
        suffix="$(echo $url | sed 's/.*\///')"

        # Fetch the genome
        wget $url/$suffix$ending
        
        gunzip $suffix$ending

        # Annotate the genome
        prokka --prefix $suffix --fast --quiet $suffix$unzip_ending
        mv $suffix/* ./

        gzip $suffix.fna
        gzip $suffix.faa

    done

    tar cvf genomes.tar *.fna.gz *.faa.gz

  }

  output {
    File genome_tar="genomes.tar"
  }

}

task PickGenomeAbundances {

  File genome_list
  String num_genomes = "100"
  String mean_depth = "1"
  String max_depth = "20"
  String sep=","
  String random_seed = "1"
  String log_std = "2"

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    set -e;

    echo "Using random seed ${random_seed}"

    python << END
import gzip
import random
import numpy as np

num_genomes = int(${num_genomes})
mean_depth = float(${mean_depth})
assert mean_depth > 0.01 and mean_depth < 100, "Please pick a reasonable depth"

# Open a file handle
f = gzip.open("${genome_list}", "rt")

# Read in the header
header = f.readline().rstrip("\n").split("${sep}")

# Figure out which column has the RefSeq accession
assert "RefSeq FTP" in header
refseq_ftp_ix = header.index("RefSeq FTP")

# Get the complete list of RefSeq FTP paths
all_genomes = [
  l.rstrip("\n").split("${sep}")[refseq_ftp_ix]
  for l in f
]
all_genomes = [g for g in all_genomes if len(g) > 1]

print("There are " + str(len(all_genomes)) + " genomes available")
assert num_genomes <= len(all_genomes)

# Randomly shuffle the genomes
random.shuffle(all_genomes)

# Randomly pick a set of abundances
log_mean = np.log10(mean_depth)
abund_list = np.random.normal(
  loc=log_mean,
  scale=${log_std},
  size=num_genomes
)

# Transform back to linear space
abund_list = [
  min(10**x, ${max_depth})
  for x in abund_list
]

# Write out to a file
with open("genome_abund.csv", "wt") as fo:
    for k, v in zip(
      all_genomes[:num_genomes],
      abund_list
    ):
        print(k, v)
        fo.write(k + ',' + str(v) + '\n')

END

  }
  output {
    File genome_abund_csv = "genome_abund.csv"
  }
}
