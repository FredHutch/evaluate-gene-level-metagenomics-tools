import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/77dcf3de35d4e2d8d11810762a37e88880001176/WDL/denovo-assembly-plass/denovo-assembly-plass.wdl" as plass
import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/3c4d7f3f5125b93981315e1c3202dce4952e5fcd/WDL/align-proteins-famli/align-proteins-famli.wdl" as famli
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/1a4791d42b50b4fb87b445a5999cb193651a6693/tools/calculate_gene_accuracy/calculate_gene_accuracy.wdl?token=AE-VSEiy1FYZRGNLtRVAC2o-hcYnRHK7ks5cXeX4wA" as calc_acc
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/1f26244639445b030282c63d60ed4f9d5365a179/tools/extract_detected_genes/extract_detected_genes_famli.wdl?token=AE-VSN7aw0nY2XpbGuT3ZEdMkFbVTqPaks5cXiLrwA" as extract_famli_genes
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/be533c85a18a7807313f84518f80a857bdfd6fd4/tools/make_unique_fasta_headers/make_unique_fasta_headers.wdl?token=AE-VSOKDLxrQh5a8n42Do3DhQNFY5tREks5cXh9LwA" as clean_headers
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/5dfcafb254884d132480a0a470d42190695c6b3c/tools/cluster_proteins_by_identity/cluster_proteins_by_identity.wdl?token=AE-VSPdWLIiGQWDCprmvobg8NMGKhksdks5cXiI8wA" as clust
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/e734a480671fe586f5efd44570bede9839c48351/tools/simulate_metagenome/simulate_metagenome.wdl?token=AE-VSHGgJIcYBtJ39Ief4u547Vo7-2ZNks5cZF3lwA" as sim

workflow evaluateGeneDetection {

  File genome_list
  Int num_metagenomes = 5
  String num_genomes = "100"
  String mean_depth = "1"
  String max_depth = "20"
  String identity="0.9"
  String memory="16G"
  String cpu="32"

  call sim.SimulateMetagenome as sim_meta {
    input:
      genome_list=genome_list,
      num_metagenomes=num_metagenomes,
      num_genomes=num_genomes,
      mean_depth=mean_depth,
      max_depth=max_depth
  }

  scatter (ix in range(num_metagenomes)){
    
    # Run Plass
    call plass.plass {
      input:
        input_fastq=sim_meta.reads_fastq[ix],
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

  output {
    Array[File] plass=plass_calc.accuracy
    Array[File] famli=famli_calc.accuracy
    Array[File] all_diamond=allDiamond_calc.accuracy
    Array[File] unique_diamond=uniqueDiamond_calc.accuracy
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