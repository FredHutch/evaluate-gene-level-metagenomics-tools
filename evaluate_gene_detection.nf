#!/usr/bin/env nextflow

genome_list_f = file(params.genome_list)
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
params.min_orf_length = 20
params.identity = 0.9
params.overlap = 0.5
params.top_pct = 0
params.output_folder = "accuracy_results/"

//
// SIMULATE MOCK COMMUNITIES
//

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
    val random_seed from params.random_seed
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
    file "reads.fastq.gz" into reads_fastq_plass

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
    file "genes.fastp.gz" into ref_genes_fasta

    script:
    template "make_gene_abundances.sh"
}

//
// RUN PLASS
//

process plass {
    // container "quay.io/biocontainers/plass@sha256:c771c791ad89d9f2c09720d7e127d5b0e6ee2a35ca7688a1b79c461c116ddd05"
    container "soedinglab/plass"
    cpus 1
    memory "8 GB"

    input:
    file input_fastq from reads_fastq_plass
    val translation_table from params.translation_table
    val min_orf_length from params.min_orf_length

    output:
    file "plass.genes.faa.gz" into plass_faa

    """
    set -e; 
    plass assemble --use-all-table-starts --min-length ${min_orf_length} --translation-table ${translation_table} "${input_fastq}" "plass.genes.faa" tmp
    gzip plass.genes.faa
    """
}

process plass_clean_headers {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 1
    memory "1 GB"

    input:
    file fasta_input from plass_faa

    output:
    file "${fasta_input}.clean.fasta.gz" into plass_clean_faa

    script:
    template "make_unique_fasta_headers.sh"
}

process plass_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 1
    memory "1 GB"

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
// CLUSTER THE REFERENCE SEQUENCES
//

process ref_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 1
    memory "1 GB"

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

process ref_cluster_abund {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 1
    memory "1 GB"

    input:
    file abund_in from gene_abund_csv
    file groups from ref_clustered_genes_tsv
    val identity from params.identity

    output:
    file "${abund_in}.clust.${identity}.csv" into ref_clustered_abund

    script:
    template "cluster_abund.sh"
}

process ref_cluster_dmnd {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 1
    memory "1 GB"

    input:
    file fasta from ref_clustered_genes_faa

    output:
    file "${fasta}.db.dmnd" into ref_clustered_genes_dmnd

    """
    diamond makedb --in ${fasta} --db ${fasta}.db.dmnd
    """
}

//
// CALCULATE THE ACCURACY OF PLASS RESULTS
//

process align_plass_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 1
    memory "1 GB"

    input:
    file db from ref_clustered_genes_dmnd
    file query from plass_clustered_faa_for_aln
    val align_id from params.identity
    val top_pct from params.top_pct
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into plass_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 1;
    gzip ${query}.${db}.aln
    """

}

process calc_plass_acc {
    container "amancevice/pandas@sha256:0c517f3aa03ac570e0cebcd2d0854f0604b44b67b7b284e79fe77307153c6f54"
    cpus 1
    memory "1 GB"
    publishDir params.output_folder

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
