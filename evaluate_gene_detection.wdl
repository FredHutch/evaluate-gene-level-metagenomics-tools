import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/0b77bb3af87598ea2e497410370032a21fca4b95/WDL/denovo-assembly-plass/denovo-assembly-plass.wdl" as plass
import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/2d3602030c21841935543c1124ecc3beca638c71/WDL/denovo-assembly-metaspades/denovo-assembly-metaspades.wdl" as metaspades
import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/2d3602030c21841935543c1124ecc3beca638c71/WDL/denovo-assembly-megahit/denovo-assembly-megahit.wdl" as megahit
import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/2d3602030c21841935543c1124ecc3beca638c71/WDL/genome-annotation-prokka/genome-annotation-prokka.wdl" as prokka
import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/3c4d7f3f5125b93981315e1c3202dce4952e5fcd/WDL/align-proteins-famli/align-proteins-famli.wdl" as famli
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/106f33a40faec5bdede8a60f3b9829f35dbe935e/tools/calculate_gene_accuracy/calculate_gene_accuracy.wdl" as calc_acc
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/106f33a40faec5bdede8a60f3b9829f35dbe935e/tools/extract_detected_genes/extract_detected_genes_famli.wdl" as extract_famli_genes
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/106f33a40faec5bdede8a60f3b9829f35dbe935e/tools/make_unique_fasta_headers/make_unique_fasta_headers.wdl" as clean_headers
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/106f33a40faec5bdede8a60f3b9829f35dbe935e/tools/cluster_proteins_by_identity/cluster_proteins_by_identity.wdl" as clust
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/106f33a40faec5bdede8a60f3b9829f35dbe935e/tools/simulate_metagenome/simulate_metagenome.wdl" as sim

workflow evaluateGeneDetection {

  File genome_list
  Int num_metagenomes = 5
  String num_genomes = "100"
  String mean_depth = "1"
  String max_depth = "20"
  String identity="0.9"
  String memory="16G"
  String cpu="32"
  String plass_min_orf_length="20"

  call sim.SimulateMetagenome as sim_meta {
    input:
      genome_list=genome_list,
      num_metagenomes=num_metagenomes,
      num_genomes=num_genomes,
      mean_depth=mean_depth,
      max_depth=max_depth
  }

  scatter (ix in range(num_metagenomes)){

    #########
    # PLASS #
    #########
    call plass.plass {
      input:
        input_fastq=sim_meta.reads_fastq[ix],
        translation_table="11",
        min_orf_length=plass_min_orf_length,
        memory=memory,
        cpu=cpu
    }
    # Make the headers unique
    call clean_headers.makeUniqueFastaHeaders as plass_clean {
      input:
        fasta_input=plass.contigs
    }
    # Cluster the Plass-detected proteins
    call clust.clusterProteins as plass_clust {
      input:
        fasta_in=plass_clean.fasta_output,
        identity=identity
    }

    # Cluster the reference proteins
    call clust.clusterProteinsAbunds as ref_clust {
      input:
        fasta_in=sim_meta.genes_fastp[ix],
        abund_in=sim_meta.abund_csv[ix],
        identity=identity
    }

    # Calculate the accuracy of Plass, with that clustering applied
    call calc_acc.calculateGeneAccuracy as plass_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=plass_clust.fasta_out,
        method_label="Plass"
    }


    ##############
    # metaSPAdes #
    ##############
    call metaspades.metaspades {
      input:
        input_fastq=sim_meta.reads_fastq[ix],
        memory=memory,
        cpu=cpu
    }
    # Annotate genes in those contigs
    call prokka.prokka as metaspades_prokka {
      input:
        input_fasta=metaspades.contigs,
        memory=memory,
        cpu=cpu
    }
    # Make the headers unique
    call clean_headers.makeUniqueFastaHeaders as metaspades_clean {
      input:
        fasta_input=metaspades_prokka.faa
    }
    # Cluster the metaSPAdes-detected proteins
    call clust.clusterProteins as metaspades_clust {
      input:
        fasta_in=metaspades_clean.fasta_output,
        identity=identity
    }
    # Calculate accuracy
    call calc_acc.calculateGeneAccuracy as metaspades_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=metaspades_clust.fasta_out,
        method_label="metaSPAdes"
    }

    ##############
    # MEGAHIT #
    ##############
    call megahit.megahit {
      input:
        input_fastq=sim_meta.reads_fastq[ix],
        memory=memory,
        cpu=cpu
    }
    # Annotate genes in those contigs
    call prokka.prokka as megahit_prokka {
      input:
        input_fasta=megahit.contigs,
        memory=memory,
        cpu=cpu
    }
    # Make the headers unique
    call clean_headers.makeUniqueFastaHeaders as megahit_clean {
      input:
        fasta_input=megahit_prokka.faa
    }
    # Cluster the megahit-detected proteins
    call clust.clusterProteins as megahit_clust {
      input:
        fasta_in=megahit_clean.fasta_output,
        identity=identity
    }
    # Calculate accuracy
    call calc_acc.calculateGeneAccuracy as megahit_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=megahit_clust.fasta_out,
        method_label="megahit"
    }

    ###########
    # DIAMOND #
    ###########
    # Now run DIAMOND on those clustered proteins
    # Make the database
    call famli.MakeDiamondDatabase {
      input:
        fasta=ref_clust.fasta_out
    }
    # Run the DIAMOND aligner
    call famli.DiamondBlastx {
      input:
        refdb=MakeDiamondDatabase.db,
        input_fastq=sim_meta.reads_fastq[ix],
        cpu=cpu
    }

    #########
    # FAMLI #
    #########
    # Filter the results with FAMLI
    call famli.FAMLI {
      input:
        input_aln=DiamondBlastx.aln,
        cpu=cpu
    }

    # Extract the genes detected by FAMLI
    call extract_famli_genes.extractDetectedGenesFAMLI as famli_genes {
      input:
        json_input=FAMLI.results,
        fasta_input=ref_clust.fasta_out
    }
    # Calculate the accuracy of the FAMLI results
    call calc_acc.calculateGeneAccuracy as famli_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=famli_genes.fasta_output,
        method_label="FAMLI"
    }

    ###############
    # All DIAMOND #
    ###############
    # Extract the set of genes which DIAMOND aligns any read to
    call allDiamondGenes {
      input:
        fasta_in=ref_clust.fasta_out,
        aln=DiamondBlastx.aln
    }
    # Calculate the accuracy of that set of genes
    call calc_acc.calculateGeneAccuracy as allDiamond_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=allDiamondGenes.fasta_output,
        method_label="All DIAMOND"
    }

    ##################
    # Unique DIAMOND #
    ##################
    # Extract the set of genes which DIAMOND aligns any read to UNIQUELY
    call uniqueDiamondGenes {
      input:
        fasta_in=ref_clust.fasta_out,
        aln=DiamondBlastx.aln
    }
    # Calculate the accuracy of that set of genes
    call calc_acc.calculateGeneAccuracy as uniqueDiamond_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=uniqueDiamondGenes.fasta_output,
        method_label="Unique DIAMOND"
    }

  }

  call AggregateResults {
    input:
      plass=plass_calc.accuracy,
      metaspades=metaspades_calc.accuracy,
      megahit=megahit_calc.accuracy,
      famli=famli_calc.accuracy,
      all_diamond=allDiamond_calc.accuracy,
      unique_diamond=uniqueDiamond_calc.accuracy
  }

  output {
    AggregateResults.results
  }

}

task AggregateResults {
  Array[File] plass
  Array[File] metaspades
  Array[File] megahit
  Array[File] famli
  Array[File] all_diamond
  Array[File] unique_diamond

  runtime {
    docker: "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    memory: "1G"
    cpu: "1"
  }

  command {
    set -e;

    python << END
import pandas as pd

plass = ['${sep="', '" plass}']
metaspades = ['${sep="', '" metaspades}']
megahit = ['${sep="', '" megahit}']
famli = ['${sep="', '" famli}']
all_diamond = ['${sep="', '" all_diamond}']
unique_diamond = ['${sep="', '" unique_diamond}']

assert len(famli) == len(plass)
assert len(all_diamond) == len(plass)
assert len(unique_diamond) == len(plass)

output = []

for ix, fp_arr in enumerate(zip(
  plass,
  metaspades,
  megahit,
  famli,
  all_diamond,
  unique_diamond
)):
    for fp in fp_arr:
        df = pd.read_table(fp)
        df["shard"] = ix
        output.append(df)

output = pd.concat(output)
output.to_csv("results.csv", index=False)
END
  }

  output {
    File results = "results.csv"
  }
}

task allDiamondGenes{
  File fasta_in
  File aln
  String sample_base = basename(aln, ".aln")

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    set -e
    python << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

detected_genes = set([])

for line in gzip_safe_open("${aln}"):
    line = line.rstrip("\n").split("\t")
    if len(line) < 5:
        continue
    detected_genes.add(line[1])

written_genes = dict()
for header, seq in SimpleFastaParser(gzip_safe_open("${fasta_in}")):
    if header in detected_genes:
        written_genes[header] = seq

assert len(written_genes) == len(detected_genes)

with open("${sample_base}.fasta", "wt") as fo:
    for k, v in written_genes.items():
        fo.write(">" + k + "\n" + v + "\n")

END
    gzip "${sample_base}.fasta"

  }

  output {
    File fasta_output = "${sample_base}.fasta.gz"
  }

}

task uniqueDiamondGenes{
  File fasta_in
  File aln
  String sample_base = basename(aln, ".aln")

  runtime {
    docker: "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    memory: "1G"
    cpu: "1"
  }

  command {
    set -e
    python << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
import gzip

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

read_counts = defaultdict(int)
read_assignments = dict()

for line in gzip_safe_open("${aln}"):
    line = line.rstrip("\n").split("\t")
    if len(line) < 5:
        continue
    read_assignments[line[0]] = line[1]
    read_counts[line[0]] += 1

detected_genes = set([
  v
  for k, v in read_assignments.items()
  if read_counts[k] == 1
])

written_genes = dict()
for header, seq in SimpleFastaParser(gzip_safe_open("${fasta_in}")):
    if header in detected_genes:
        written_genes[header] = seq

assert len(written_genes) == len(detected_genes)

with open("${sample_base}.fasta", "wt") as fo:
    for k, v in written_genes.items():
        fo.write(">" + k + "\n" + v + "\n")

END
    gzip "${sample_base}.fasta"

  }

  output {
    File fasta_output = "${sample_base}.fasta.gz"
  }

}