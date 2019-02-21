#!/usr/bin/env nextflow

genome_list_f = file(params.genome_list)
params.genome_list_sep = ","
params.num_metagenomes = 10
params.num_genomes = 10
params.mean_depth = 5
params.max_depth = 20
params.log_std = 2
params.random_seed = 1

metagenome_seed_ch = Channel.from( 0..params.num_metagenomes)

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
    file "genome_abund.csv" into genome_abund_csv

    script:
    template "pick_genome_abundances.sh"

}

process download_genomes {

    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 1
    memory "4 GB"

    input:
    file genome_abund_csv

    output:
    file "genomes.tar"

    script:
    template "download_genomes.sh"

}
