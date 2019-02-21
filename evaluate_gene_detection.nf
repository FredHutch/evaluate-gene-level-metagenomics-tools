#!/usr/bin/env nextflow

genome_list_f = file(params.genome_list)
params.genome_list_sep = ","
params.num_metagenomes = 2
params.num_genomes = 2
params.mean_depth = 5
params.max_depth = 2
params.log_std = 2
params.read_length = 100
params.read_mflen = 300
params.read_sdev = 75

metagenome_seed_ch = Channel.from( 1..params.num_metagenomes)

process pick_genome_abundances {

    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 1
    memory "1 GB"

    input:
    file genome_list from genome_list_f
    val num_genomes from params.num_genomes
    val mean_depth from params.mean_depth
    val max_depth from params.max_depth
    val log_std from params.log_std
    val random_seed from metagenome_seed_ch
    val sep from params.genome_list_sep

    output:
    file "genome_abund.csv" into genome_abund_csv_dg

    script:
    template "pick_genome_abundances.sh"

}

process download_genomes {

    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 1
    memory "1 GB"

    input:
    file genome_abund_csv from genome_abund_csv_dg

    output:
    file "genomes.tar" into genome_tar_sg, genome_tar_mga
    file "genome_abund.csv" into genome_abund_csv_sg, genome_abund_csv_mga

    script:
    template "download_genomes.sh"

}

process simulate_genomes {
    container "quay.io/biocontainers/art@sha256:14f44c1cf099f6b55922aaa177c926c993733dba75ef4eb4dcf53442e3b5f96e"
    cpus 1
    memory "1 GB"

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

process interleave_fastqs {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 1
    memory "1 GB"

    input:
    file all_reads_tar

    output:
    file "reads.fastq.gz" into reads_fastq

    script:
    template "interleave_fastq.sh"

}

process make_gene_abundances {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 1
    memory "1 GB"

    input:
    file genome_abund_csv from genome_abund_csv_mga
    file genome_tar from genome_tar_mga

    output:
    file "genes_abund.csv.gz" into gene_abund_csv
    file "genes.fastp.gz" into genes_fastp

    script:
    template "make_gene_abundances.sh"
}