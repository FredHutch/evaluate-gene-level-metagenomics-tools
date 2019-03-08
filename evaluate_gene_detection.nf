#!/usr/bin/env nextflow

genome_list_f = file(params.genome_list)
refdb_ch = file(params.refdb)
params.random_seed = 1
params.genome_list_sep = ","
params.num_genomes = 20
params.mean_depth = 5
params.max_depth = 50
params.log_std = 2
params.read_length = 100
params.read_mflen = 300
params.read_sdev = 75
params.translation_table = 11
params.min_orf_length = 30
params.identity = 0.9
params.overlap = 0.5
params.output_folder = "accuracy_results/"

//
// SIMULATE MOCK COMMUNITIES
//

process pick_genome_abundances {

    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file genome_list from genome_list_f
    val num_genomes from params.num_genomes
    val mean_depth from params.mean_depth
    val max_depth from params.max_depth
    val log_std from params.log_std
    val random_seed from params.random_seed
    val sep from params.genome_list_sep

    output:
    file "genome_abund.${random_seed}.csv" into genome_abund_csv_dg

    script:
    template "pick_genome_abundances.sh"

}

//
// DOWNLOAD GENOMES FOR SIMULATION
//

process download_genomes {

    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 16
    memory "32 GB"
    scratch "/scratch"
    publishDir params.output_folder

    input:
    file genome_abund_csv from genome_abund_csv_dg
    val random_seed from params.random_seed

    output:
    file "genomes.${random_seed}.tar" into genome_tar_sg, genome_tar_mga
    file "genome_abund.${random_seed}.csv" into genome_abund_csv_sg, genome_abund_csv_mga

    script:
    template "download_genomes.sh"

}

//
// SIMULATE READS FROM DOWNLOADED GENOMES
//

process simulate_genomes {
    container "quay.io/biocontainers/art@sha256:14f44c1cf099f6b55922aaa177c926c993733dba75ef4eb4dcf53442e3b5f96e"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    val read_length from params.read_length
    val mflen from params.read_mflen
    val sdev from params.read_sdev
    file genome_abund_csv from genome_abund_csv_sg
    file genome_tar from genome_tar_sg

    output:
    file "all_reads.tar" into all_reads_tar

    script:
    template "simulate_genomes.sh"
}

//
// INTERLEAVE THE READS FROM EACH SIMULATED GENOME
//

process interleave_fastqs {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file all_reads_tar

    output:
    file "reads.fastq.gz" into reads_fastq_plass, reads_fastq_metaspades, reads_fastq_megahit, reads_fastq_diamond

    script:
    template "interleave_fastq.sh"

}

//
// CALCULATE THE ABUNDANCES OF EACH INDIVIDUAL GENE IN THE SIMULATION
//

process make_gene_abundances {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    scratch "/scratch"
    publishDir params.output_folder

    input:
    file genome_abund_csv from genome_abund_csv_mga
    file genome_tar from genome_tar_mga
    val random_seed from params.random_seed

    output:
    file "genes_abund.${random_seed}.csv.gz" into gene_abund_csv
    file "genes.${random_seed}.fastp.gz" into ref_genes_fasta

    script:
    template "make_gene_abundances.sh"
}

//
// RUN PLASS
//

process plass {
    container "quay.io/fhcrc-microbiome/plass@sha256:72d2c563a7ed97c20064116656f93edbb7c92d0cce7ee4a9f5189fcbbbcad13f"
    cpus 16
    memory "122 GB"
    scratch "/scratch"

    input:
    file input_fastq from reads_fastq_plass
    val translation_table from params.translation_table
    val min_orf_length from params.min_orf_length

    output:
    file "plass.genes.faa.gz" into plass_faa

    """
    set -e; 
    /usr/local/plass/build/bin/plass assemble --use-all-table-starts --min-length ${min_orf_length} --threads 16 --translation-table ${translation_table} "${input_fastq}" "plass.genes.faa" tmp
    gzip plass.genes.faa
    """
}

//
// CLEAN UP THE HEADERS OF THE PLASS OUTPUT
//

process plass_clean_headers {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    cpus 1
    memory "4 GB"
    scratch "/scratch"

    input:
    file fasta_input from plass_faa

    output:
    file "${fasta_input}.clean.fasta.gz" into plass_clean_faa

    script:
    template "make_unique_fasta_headers.sh"
}

//
// CLUSTER THE PLASS GENES
//

process plass_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file fasta_in from plass_clean_faa
    val identity from params.identity
    val overlap from params.overlap

    output:
    file "${fasta_in}.rep.fasta" into plass_clustered_faa_for_aln, plass_clustered_faa_for_acc
    file "${fasta_in}.clusters.tsv" into plass_clustered_tsv

    script:
    template "cluster_proteins.sh"

}

//
// CLUSTER THE REFERENCE GENES
//

process ref_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file fasta_in from ref_genes_fasta
    val identity from params.identity
    val overlap from params.overlap

    output:
    file "${fasta_in}.rep.fasta" into ref_clustered_genes_faa
    file "${fasta_in}.clusters.tsv" into ref_clustered_genes_tsv

    script:
    template "cluster_proteins.sh"

}

//
// CALCULATE THE ABUNDANCES OF THE CLUSTERED REFERENCE GENES
//

process ref_cluster_abund {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file abund_in from gene_abund_csv
    file groups from ref_clustered_genes_tsv
    val identity from params.identity

    output:
    file "${abund_in}.clust.${identity}.csv" into ref_clustered_abund

    script:
    template "cluster_abund.sh"
}

//
// MAKE A DIAMOND DATABASE FOR THE CLUSTERED REFERENCE GENES
//

process ref_cluster_dmnd {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file fasta from ref_clustered_genes_faa

    output:
    file "${fasta}.db.dmnd" into ref_clustered_genes_dmnd

    """
    diamond makedb --in ${fasta} --db ${fasta}.db.dmnd
    """
}

//
// ALIGN THE PLASS GENES AGAINST THE REFERENCE GENES
//

process align_plass_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file db from ref_clustered_genes_dmnd
    file query from plass_clustered_faa_for_aln
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into plass_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE THE ACCURACY OF PLASS RESULTS
//

process calc_plass_acc {
    container "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    cpus 16
    memory "32 GB"
    publishDir params.output_folder
    scratch "/scratch"

    input:
    file detected_fasta from plass_clustered_faa_for_acc
    file aln from plass_ref_aln
    file ref_abund from ref_clustered_abund
    val method_label from "Plass"
    val random_seed from params.random_seed

    output:
    file "${method_label}.${random_seed}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// RUN METASPADES
//

process metaspades {
    container "quay.io/biocontainers/spades@sha256:9f097c5d6d7944b68828e10d94504ac49a93bf337a9afed17232594b126b807e"
    cpus 16
    memory "122 GB"
    scratch "/scratch"

    input:
    file input_fastq from reads_fastq_metaspades

    output:
    file "${input_fastq}.metaspades.fasta.gz" into metaspades_contigs_fasta

    """
    set -e; 

    spades.py --12 ${input_fastq} -o TEMP --threads 16 --meta --phred-offset 33 --only-assembler;
    
    mv TEMP/scaffolds.fasta ${input_fastq}.metaspades.fasta
    gzip ${input_fastq}.metaspades.fasta
    """
}

//
// GET METASPADES GENES
//

process metaspades_prokka {
    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file input_fasta from metaspades_contigs_fasta

    output:
    file "TEMP/metaspades.faa.gz" into metaspades_genes_fasta

    """
    set -e; 
    
    gunzip -c ${input_fasta} > input.fasta; 
    prokka --outdir TEMP --prefix metaspades --metagenome input.fasta
    gzip TEMP/metaspades.faa
    """

}

//
// CLUSTER METASPADES GENES
//

process metaspades_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file fasta_in from metaspades_genes_fasta
    val identity from params.identity
    val overlap from params.overlap

    output:
    file "${fasta_in}.rep.fasta" into metaspades_clustered_faa_for_aln, metaspades_clustered_faa_for_acc
    file "${fasta_in}.clusters.tsv" into metaspades_clustered_tsv

    script:
    template "cluster_proteins.sh"

}

//
// ALIGN METASPADES GENES AGAINST THE REFERENCE GENES
//


process align_metaspades_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file db from ref_clustered_genes_dmnd
    file query from metaspades_clustered_faa_for_aln
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into metaspades_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE ACCURACY OF METASPADES GENES
//

process calc_metaspades_acc {
    container "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    cpus 16
    memory "32 GB"
    publishDir params.output_folder
    scratch "/scratch"

    input:
    file detected_fasta from metaspades_clustered_faa_for_acc
    file aln from metaspades_ref_aln
    file ref_abund from ref_clustered_abund
    val method_label from "metaSPAdes"
    val random_seed from params.random_seed

    output:
    file "${method_label}.${random_seed}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// RUN MEGAHIT
//

process megahit {
    container "quay.io/biocontainers/megahit@sha256:8c9f17dd0fb144254e4d6a2a11d46b522239d752d2bd15ae3053bb1a31cc6d01"
    cpus 16
    memory "122 GB"
    scratch "/scratch"

    input:
    file input_fastq from reads_fastq_megahit

    output:
    file "megahit.contigs.fasta.gz" into megahit_contigs_fasta

    """
    set -e
    megahit --12 ${input_fastq} -o TEMP -t 16
    mv TEMP/final.contigs.fa megahit.contigs.fasta
    [[ -s megahit.contigs.fasta ]]
    gzip megahit.contigs.fasta
    """
}

//
// GET MEGAHIT GENES
//

process megahit_prokka {
    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file input_fasta from megahit_contigs_fasta

    output:
    file "TEMP/megahit.faa.gz" into megahit_genes_fasta

    """
    set -e; 
    
    gunzip -c ${input_fasta} > input.fasta; 
    prokka --outdir TEMP --prefix megahit --metagenome input.fasta
    gzip TEMP/megahit.faa
    """

}

//
// CLUSTER MEGAHIT GENES
//

process megahit_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file fasta_in from megahit_genes_fasta
    val identity from params.identity
    val overlap from params.overlap

    output:
    file "${fasta_in}.rep.fasta" into megahit_clustered_faa_for_aln, megahit_clustered_faa_for_acc
    file "${fasta_in}.clusters.tsv" into megahit_clustered_tsv

    script:
    template "cluster_proteins.sh"

}

//
// ALIGN MEGAHIT GENES AGAINST THE REFERENCE GENES
//

process align_megahit_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file db from ref_clustered_genes_dmnd
    file query from megahit_clustered_faa_for_aln
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into megahit_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE ACCURACY OF MEGAHIT GENES
//

process calc_megahit_acc {
    container "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    cpus 16
    memory "32 GB"
    publishDir params.output_folder
    scratch "/scratch"

    input:
    file detected_fasta from megahit_clustered_faa_for_acc
    file aln from megahit_ref_aln
    file ref_abund from ref_clustered_abund
    val method_label from "megahit"
    val random_seed from params.random_seed

    output:
    file "${method_label}.${random_seed}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// RUN DIAMOND TO ALIGN READS AGAINST REFERENCE PROTEINS
//

process make_refdb_dmnd {
    // container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    container "job-definition://diamond_nf:3"
    cpus 16
    memory "120 GB"
    scratch "/scratch"
    
    input:
    file fasta from refdb_ch

    output:
    file "${fasta}.db.dmnd" into refdb_dmnd

    """
    diamond makedb --in ${fasta} --db ${fasta}.db.dmnd
    """
}

process diamond {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "120 GB"
    scratch "/scratch"

    input:
    file refdb from refdb_dmnd
    file input_fastq from reads_fastq_diamond
    val min_id from params.identity
    val query_cover from params.overlap
    val subject_cover from params.overlap
    val cpu from 4
    val min_score from 20
    val blocks from 20
    val query_gencode from 11

    output:
    file "${input_fastq}.${refdb}.aln" into diamond_aln

    """
    set -e
    diamond \
      blastx \
      --query ${input_fastq} \
      --out ${input_fastq}.${refdb}.aln \
      --threads ${cpu} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${min_score} \
      --query-cover ${query_cover * 100} \
      --id ${min_id * 100} \
      --top ${100 - (100 * min_id)} \
      --block-size ${blocks} \
      --query-gencode ${query_gencode} \
      --unal 0    
    """

}

//
// RUN FAMLI TO FILTER DIAMOND ALIGNMENTS
//

process famli {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file input_aln from diamond_aln
    val cpu from 1
    val batchsize from 50000000

    output:
    file "${input_aln}.json.gz" into famli_json

    """
    set -e; 
    
    famli \
      filter \
      --input ${input_aln} \
      --output ${input_aln}.json \
      --threads ${cpu} \
      --batchsize ${batchsize}

    gzip ${input_aln}.json
    """
}

//
// GET GENES DETECTED BY FAMLI
//

process famli_genes {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file json_input from famli_json
    file fasta_input from refdb_ch

    output:
    file "${json_input}.faa.gz" into famli_clustered_faa_for_acc, famli_clustered_faa_for_aln

    script:
    template "extract_detected_genes_famli.sh"
}

//
// ALIGN FAMLI GENES AGAINST THE REFERENCE GENES
//

process align_famli_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "120 GB"
    scratch "/scratch"

    input:
    file db from ref_clustered_genes_dmnd
    file query from famli_clustered_faa_for_aln
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into famli_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE THE ACCURACY OF THE FAMLI RESULTS
//

process calc_famli_acc {
    container "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    cpus 16
    memory "32 GB"
    publishDir params.output_folder
    scratch "/scratch"

    input:
    file detected_fasta from famli_clustered_faa_for_acc
    file aln from famli_ref_aln
    file ref_abund from ref_clustered_abund
    val method_label from "FAMLI"
    val random_seed from params.random_seed

    output:
    file "${method_label}.${random_seed}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// GET THE TOTAL SET OF GENES WITH ANY READ ALIGNING BY DIAMOND
//

process all_diamond_genes {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file fasta_in from refdb_ch
    file aln from diamond_aln

    output:
    file "${aln}.fasta.gz" into all_diamond_faa_for_acc, all_diamond_faa_for_aln

    script:
    template "all_diamond_genes.sh"
}

//
// ALIGN ALL-DIAMOND GENES AGAINST THE REFERENCE GENES
//

process align_all_diamond_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file db from ref_clustered_genes_dmnd
    file query from all_diamond_faa_for_aln
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into all_diamond_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE THE ACCURACY OF THE ALL-DIAMOND RESULTS
//

process calc_all_diamond_acc {
    container "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    cpus 16
    memory "32 GB"
    publishDir params.output_folder
    scratch "/scratch"

    input:
    file detected_fasta from all_diamond_faa_for_acc
    file aln from all_diamond_ref_aln
    file ref_abund from ref_clustered_abund
    val method_label from "all_diamond"
    val random_seed from params.random_seed

    output:
    file "${method_label}.${random_seed}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// GET THE TOTAL SET OF GENES WITH ANY READ ALIGNING *UNIQUELY* BY DIAMOND
//

process unique_diamond_genes {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file fasta_in from refdb_ch
    file aln from diamond_aln

    output:
    file "${aln}.fasta.gz" into unique_diamond_faa_for_acc, unique_diamond_faa_for_aln

    script:
    template "unique_diamond_genes.sh"
}

//
// ALIGN UNIQUE-DIAMOND GENES AGAINST THE REFERENCE GENES
//

process align_unique_diamond_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "32 GB"
    scratch "/scratch"

    input:
    file db from ref_clustered_genes_dmnd
    file query from unique_diamond_faa_for_aln
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into unique_diamond_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE THE ACCURACY OF THE UNIQUE-DIAMOND RESULTS
//

process calc_unique_diamond_acc {
    container "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    cpus 16
    memory "32 GB"
    publishDir params.output_folder
    scratch "/scratch"

    input:
    file detected_fasta from unique_diamond_faa_for_acc
    file aln from unique_diamond_ref_aln
    file ref_abund from ref_clustered_abund
    val method_label from "unique_diamond"
    val random_seed from params.random_seed

    output:
    file "${method_label}.${random_seed}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}