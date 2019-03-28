#!/usr/bin/env nextflow

genome_list_f = file(params.genome_list)
refdb_f = file(params.refdb)
params.num_simulations = 100
params.genome_list_sep = ","
params.num_genomes = 20
params.mean_depth = 5
params.max_depth = 100
params.log_std = 1
params.read_length = 250
params.read_mflen = 1000
params.read_sdev = 300
params.translation_table = 11
params.min_orf_length = 30
params.identity = 0.9
params.overlap = 0.5
params.output_folder = "accuracy_results/"

// Make a channel with the index for each simulation
simulation_index_ch = Channel.from( 1..params.num_simulations )

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
    val ix from simulation_index_ch
    val sep from params.genome_list_sep

    output:
    set ix, "genome_abund.${ix}.csv" into genome_abund_csv_dg

    script:
    template "pick_genome_abundances.sh"

}

//
// DOWNLOAD GENOMES FOR SIMULATION
//

process download_genomes {

    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 4
    memory "8 GB"
    publishDir params.output_folder

    input:
    set ix, "genome_abund.${ix}.csv" from genome_abund_csv_dg

    output:
    set ix, "genome_abund.${ix}.csv", "genomes.${ix}.tar" into genome_tar_sg, genome_tar_mga

    script:
    template "download_genomes.sh"

}

//
// SIMULATE READS FROM DOWNLOADED GENOMES
//

process simulate_genomes {
    container "quay.io/biocontainers/art@sha256:1cd93ed9f680318812a9b19b6187e24c8a63e14cd0fd0202844bc46a8d9ed2b2"
    cpus 4
    memory "8 GB"
    errorStrategy 'retry'

    input:
    val read_length from params.read_length
    val mflen from params.read_mflen
    val sdev from params.read_sdev
    set ix, "genome_abund.${ix}.csv", "genomes.${ix}.tar" from genome_tar_sg

    output:
    set ix, "all_reads.${ix}.tar" into all_reads_tar

    script:
    template "simulate_genomes.sh"
}

//
// INTERLEAVE THE READS FROM EACH SIMULATED GENOME
//

process interleave_fastqs {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 4
    memory "8 GB"

    input:
    set ix, "all_reads.${ix}.tar" from all_reads_tar

    output:
    set ix, "reads.${ix}.fastq.gz" into reads_fastq_plass, reads_fastq_metaspades, reads_fastq_megahit, reads_fastq_idba, reads_fastq_diamond, reads_fastq_humann2

    script:
    template "interleave_fastq.sh"

}

//
// CALCULATE THE ABUNDANCES OF EACH INDIVIDUAL GENE IN THE SIMULATION
//

process make_gene_abundances {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 1
    memory "1 GB"
    publishDir params.output_folder

    input:
    set ix, "genome_abund.${ix}.csv", "genomes.${ix}.tar" from genome_tar_mga
    
    output:
    set ix, "genes_abund.${ix}.csv.gz" into gene_abund_csv
    set ix, "genes.${ix}.fastp.gz" into ref_genes_fasta

    script:
    template "make_gene_abundances.sh"
}

//
// RUN PLASS
//

process plass {
    container "quay.io/fhcrc-microbiome/plass@sha256:72d2c563a7ed97c20064116656f93edbb7c92d0cce7ee4a9f5189fcbbbcad13f"
    cpus 32
    memory "240 GB"

    input:
    set ix, file(fastq_in) from reads_fastq_plass
    val translation_table from params.translation_table
    val min_orf_length from params.min_orf_length

    output:
    set ix, "plass.${ix}.genes.faa.gz" into plass_faa

    """
    set -e; 
    /usr/local/plass/build/bin/plass assemble --use-all-table-starts --min-length ${min_orf_length} --threads 16 --translation-table ${translation_table} "${fastq_in}" "plass.${ix}.genes.faa" tmp
    gzip plass.${ix}.genes.faa
    """
}

//
// CLEAN UP THE HEADERS OF THE PLASS OUTPUT
//

process plass_clean_headers {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 1
    memory "4 GB"

    input:
    set ix, file(fasta_input) from plass_faa
    
    output:
    set ix, "${fasta_input}.clean.fasta.gz" into plass_clean_faa

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

    input:
    set ix, file(fasta_in) from plass_clean_faa
    val identity from params.identity
    val overlap from params.overlap

    output:
    set ix, "${fasta_in}.rep.fasta" into plass_clustered_faa_for_aln, plass_clustered_faa_for_acc
    set ix, "${fasta_in}.clusters.tsv" into plass_clustered_tsv

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

    input:
    set ix, file(fasta_in) from ref_genes_fasta
    val identity from params.identity
    val overlap from params.overlap

    output:
    set ix, "${fasta_in}.rep.fasta" into ref_clustered_genes_faa
    set ix, "${fasta_in}.clusters.tsv" into ref_clustered_genes_tsv

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

    input:
    set ix, file(abund_in), file(groups) from gene_abund_csv.join(ref_clustered_genes_tsv)
    val identity from params.identity

    output:
    set ix, "${abund_in}.clust.${identity}.csv" into ref_clustered_abund_metaspades, ref_clustered_abund_megahit, ref_clustered_abund_idba, ref_clustered_abund_plass, ref_clustered_abund_unique_diamond, ref_clustered_abund_all_diamond, ref_clustered_abund_famli, ref_clustered_abund_humann2

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

    input:
    set ix, file(fasta) from ref_clustered_genes_faa

    output:
    set ix, "${fasta}.db.dmnd" into ref_clustered_genes_dmnd_plass, ref_clustered_genes_dmnd_metaspades, ref_clustered_genes_dmnd_megahit, ref_clustered_genes_dmnd_idba, ref_clustered_genes_dmnd_famli, ref_clustered_genes_dmnd_all_diamond, ref_clustered_genes_dmnd_unique_diamond, ref_clustered_genes_dmnd_humann2

    """
    diamond makedb --in ${fasta} --db ${fasta}.db.dmnd
    """
}

//
// ALIGN THE PLASS GENES AGAINST THE REFERENCE GENES
//

process align_plass_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"

    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_plass.join(plass_clustered_faa_for_aln)
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    set ix, "${query}.${db}.aln.gz" into plass_ref_aln

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
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder

    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from plass_clustered_faa_for_acc.join(plass_ref_aln).join(ref_clustered_abund_plass)
    val method_label from "Plass"

    output:
    file "${method_label}.${ix}.accuracy.tsv"

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

    input:
    set ix, file(input_fastq) from reads_fastq_metaspades

    output:
    set ix, "${input_fastq}.metaspades.fasta.gz" into metaspades_contigs_fasta

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

    input:
    set ix, file(input_fasta) from metaspades_contigs_fasta

    output:
    set ix, "TEMP/metaspades.${ix}.faa.gz" into metaspades_genes_fasta

    """
    set -e; 
    
    gunzip -c ${input_fasta} > input.fasta; 
    prokka --outdir TEMP --prefix metaspades.${ix} --metagenome input.fasta
    gzip TEMP/metaspades.${ix}.faa
    """

}

//
// CLUSTER METASPADES GENES
//

process metaspades_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 16
    memory "32 GB"

    input:
    set ix, file(fasta_in) from metaspades_genes_fasta
    val identity from params.identity
    val overlap from params.overlap

    output:
    set ix, "${fasta_in}.rep.fasta" into metaspades_clustered_faa_for_aln, metaspades_clustered_faa_for_acc
    set ix, "${fasta_in}.clusters.tsv" into metaspades_clustered_tsv

    script:
    template "cluster_proteins.sh"

}

//
// ALIGN METASPADES GENES AGAINST THE REFERENCE GENES
//


process align_metaspades_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"

    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_metaspades.join(metaspades_clustered_faa_for_aln)
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    set ix, "${query}.${db}.aln.gz" into metaspades_ref_aln

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
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder

    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from metaspades_clustered_faa_for_acc.join(metaspades_ref_aln).join(ref_clustered_abund_metaspades)
    val method_label from "metaSPAdes"

    output:
    set ix, "${method_label}.${ix}.accuracy.tsv"

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

    input:
    set ix, file(input_fastq) from reads_fastq_megahit

    output:
    set ix, "megahit.${ix}.contigs.fasta.gz" into megahit_contigs_fasta

    """
    set -e
    megahit --12 ${input_fastq} -o TEMP -t 16
    mv TEMP/final.contigs.fa megahit.${ix}.contigs.fasta
    [[ -s megahit.${ix}.contigs.fasta ]]
    gzip megahit.${ix}.contigs.fasta
    """
}

//
// GET MEGAHIT GENES
//

process megahit_prokka {
    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 16
    memory "32 GB"

    input:
    set ix, file(input_fasta) from megahit_contigs_fasta

    output:
    set ix, "TEMP/megahit.${ix}.faa.gz" into megahit_genes_fasta

    """
    set -e; 
    
    gunzip -c ${input_fasta} > input.fasta; 
    prokka --outdir TEMP --prefix megahit.${ix} --metagenome input.fasta
    gzip TEMP/megahit.${ix}.faa
    """

}

//
// CLUSTER MEGAHIT GENES
//

process megahit_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 16
    memory "32 GB"

    input:
    set ix, file(fasta_in) from megahit_genes_fasta
    val identity from params.identity
    val overlap from params.overlap

    output:
    set ix, "${fasta_in}.rep.fasta" into megahit_clustered_faa_for_aln, megahit_clustered_faa_for_acc
    set ix, "${fasta_in}.clusters.tsv" into megahit_clustered_tsv

    script:
    template "cluster_proteins.sh"

}

//
// ALIGN MEGAHIT GENES AGAINST THE REFERENCE GENES
//

process align_megahit_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"

    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_megahit.join(megahit_clustered_faa_for_aln)
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
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder

    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from megahit_clustered_faa_for_acc.join(megahit_ref_aln).join(ref_clustered_abund_megahit)
    val method_label from "megahit"

    output:
    file "${method_label}.${ix}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// RUN IDBA-UD
//

process idba {
    container "quay.io/biocontainers/idba@sha256:51291ffeeecc6afab8d56bf33dffd0c2cb5e24d8a545a5ea93bb795d6af12fa0"
    cpus 16
    memory "122 GB"

    input:
    set ix, file(input_fastq) from reads_fastq_idba

    output:
    set ix, "idba.${ix}.scaffold.fasta.gz" into idba_contigs_fasta

    """
    set -e
    fq2fa --paired --filter <(gunzip -c ${input_fastq}) ${input_fastq}.fa
    idba_ud --num_threads 16 -r ${input_fastq}.fa -o TEMP
    mv TEMP/scaffold.fa idba.${ix}.scaffold.fasta
    [[ -s idba.${ix}.scaffold.fasta ]]
    gzip idba.${ix}.scaffold.fasta
    """
}

//
// GET IDBA GENES
//

process idba_prokka {
    container "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    cpus 16
    memory "32 GB"

    input:
    set ix, file(input_fasta) from idba_contigs_fasta

    output:
    set ix, "TEMP/idba.${ix}.faa.gz" into idba_genes_fasta

    """
    set -e; 
    
    gunzip -c ${input_fasta} > input.fasta; 
    prokka --outdir TEMP --prefix idba.${ix} --metagenome input.fasta
    gzip TEMP/idba.${ix}.faa
    """

}

//
// CLUSTER IDBA GENES
//

process idba_cluster {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    cpus 16
    memory "32 GB"

    input:
    set ix, file(fasta_in) from idba_genes_fasta
    val identity from params.identity
    val overlap from params.overlap

    output:
    set ix, "${fasta_in}.rep.fasta" into idba_clustered_faa_for_aln, idba_clustered_faa_for_acc
    set ix, "${fasta_in}.clusters.tsv" into idba_clustered_tsv

    script:
    template "cluster_proteins.sh"

}

//
// ALIGN IDBA GENES AGAINST THE REFERENCE GENES
//

process align_idba_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"

    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_idba.join(idba_clustered_faa_for_aln)
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    file "${query}.${db}.aln.gz" into idba_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE ACCURACY OF IDBA GENES
//

process calc_idba_acc {
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder

    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from idba_clustered_faa_for_acc.join(idba_ref_aln).join(ref_clustered_abund_idba)
    val method_label from "IDBA"

    output:
    set ix, "${method_label}.${ix}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// RUN DIAMOND TO ALIGN READS AGAINST REFERENCE PROTEINS
//

process make_refdb_dmnd {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 16
    memory "120 GB"
    
    input:
    file fasta from refdb_f

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

    input:
    file refdb from refdb_dmnd
    set ix, file(input_fastq) from reads_fastq_diamond
    val min_id from params.identity
    val query_cover from params.overlap
    val subject_cover from params.overlap
    val cpu from 4
    val min_score from 20
    val blocks from 20
    val query_gencode from 11

    output:
    set ix, "${input_fastq}.${refdb}.aln" into diamond_aln_all, diamond_aln_unique, diamond_aln_famli

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
    memory "120 GB"

    input:
    set ix, file(input_aln) from diamond_aln_famli
    val cpu from 1
    val batchsize from 50000000

    output:
    set ix, "${input_aln}.json.gz" into famli_json

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

    input:
    set ix, file(json_input) from famli_json
    file fasta_input from refdb_f

    output:
    set ix, "${json_input}.faa.gz" into famli_clustered_faa_for_acc, famli_clustered_faa_for_aln

    script:
    template "extract_detected_genes_famli.sh"
}

//
// ALIGN FAMLI GENES AGAINST THE REFERENCE GENES
//

process align_famli_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"

    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_famli.join(famli_clustered_faa_for_aln)
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    set ix, "${query}.${db}.aln.gz" into famli_ref_aln

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
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder

    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from famli_clustered_faa_for_acc.join(famli_ref_aln).join(ref_clustered_abund_famli)
    val method_label from "FAMLI"

    output:
    file "${method_label}.${ix}.accuracy.tsv"

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

    input:
    file fasta_in from refdb_f
    set ix, file(aln) from diamond_aln_all

    output:
    set ix, "${aln}.fasta.gz" into all_diamond_faa_for_acc, all_diamond_faa_for_aln

    script:
    template "all_diamond_genes.sh"
}

//
// ALIGN ALL-DIAMOND GENES AGAINST THE REFERENCE GENES
//

process align_all_diamond_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"

    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_all_diamond.join(all_diamond_faa_for_aln)
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    set ix, "${query}.${db}.aln.gz" into all_diamond_ref_aln

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
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder

    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from all_diamond_faa_for_acc.join(all_diamond_ref_aln).join(ref_clustered_abund_all_diamond)
    val method_label from "all_diamond"

    output:
    file "${method_label}.${ix}.accuracy.tsv"

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

    input:
    file fasta_in from refdb_f
    set ix, file(aln) from diamond_aln_unique

    output:
    set ix, "${aln}.fasta.gz" into unique_diamond_faa_for_acc, unique_diamond_faa_for_aln

    script:
    template "unique_diamond_genes.sh"
}

//
// ALIGN UNIQUE-DIAMOND GENES AGAINST THE REFERENCE GENES
//

process align_unique_diamond_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"

    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_unique_diamond.join(unique_diamond_faa_for_aln)
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
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder

    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from unique_diamond_faa_for_acc.join(unique_diamond_ref_aln).join(ref_clustered_abund_unique_diamond)
    val method_label from "unique_diamond"

    output:
    file "${method_label}.${ix}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}

//
// MAKE HUMANN2 DATABASE

process humann2_db {

    container "quay.io/fhcrc-microbiome/humann2@sha256:d6426bda36ca6a689ea7ddc1fd8c628c6e036d90469234ac1379fa9a4f4d1840"
    cpus 16
    memory "120 GB"

    output:
    file "chocophlan/*" into humann2_chocophlan
    file "uniref90_ec_filtered_diamond/*" into humann2_uniref90
    
    """
    # Download the database
	humann2_databases --download chocophlan full chocophlan
	humann2_databases --download uniref uniref90_ec_filtered_diamond uniref90_ec_filtered_diamond
	"""
}

//
// RUN HUMANN2
//

process humann2 {
    container "quay.io/fhcrc-microbiome/humann2@sha256:d6426bda36ca6a689ea7ddc1fd8c628c6e036d90469234ac1379fa9a4f4d1840"
    cpus 16
    memory "120 GB"

    input:
    set ix, file(fastq) from reads_fastq_humann2
    file chocophlan from humann2_chocophlan
    file uniref90 from humann2_uniref90
    
    output:
    set ix, "OUTPUT_DIR/${fastq.simpleName}.${ix}_genefamilies.tsv" into humann2_gene_tsv

    """
    humann2 \
    --input ${fastq} \
    --nucleotide-database ${chocophlan} \
    --protein-database ${uniref90} \
    --output OUTPUT_DIR

    rm -r ${chocophlan} ${uniref90}
    """
}

//
// GET HUMANN2 GENES
//

process extract_humann2_db {
    container "quay.io/fhcrc-microbiome/humann2@sha256:d6426bda36ca6a689ea7ddc1fd8c628c6e036d90469234ac1379fa9a4f4d1840"
    cpus 16
    memory "32 GB"

    input:
    file uniref90 from humann2_uniref90

    output:
    file "humann2_refdb.fasta" into humann2_refdb_fasta

    """
    diamond getseq -d \$(find . -name "*dmnd") > humann2_refdb.fasta
    """
}

process extract_humann2_genes {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    
    input:
    file ref_fasta from humann2_refdb_fasta
    set ix, file(gene_tsv) from humann2_gene_tsv
    file chocophlan from humann2_chocophlan

    output:
    set ix, "humann2.${ix}.detected.fasta" into humann_faa_for_aln, humann_faa_for_acc

    script:
    template "extract_detected_genes_humann2.sh"
}

//
// ALIGN HUMANN2-DETECTED GENES AGAINST THE REFERENCE GENES
//

process align_humann2_ref {
    container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    cpus 4
    memory "8 GB"
    
    input:
    set ix, file(db), file(query) from ref_clustered_genes_dmnd_humann2.join(humann_faa_for_aln)
    val align_id from params.identity
    val top_pct from 0
    val query_cover from params.overlap
    val subject_cover from params.overlap

    output:
    set ix, "${query}.${db}.aln.gz" into humann2_ref_aln

    """
    set -e;
    diamond blastp --db ${db} --query ${query} --out ${query}.${db}.aln --outfmt 6 --id ${align_id * 100} --top ${top_pct} --query-cover ${query_cover * 100} --subject-cover ${subject_cover * 100} --threads 16;
    gzip ${query}.${db}.aln
    """

}

//
// CALCULATE THE ACCURACY OF THE HUMANN2 RESULTS
//

process calc_humann2_acc {
    container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
    cpus 1
    memory "2 GB"
    publishDir params.output_folder
 
    input:
    set ix, file(detected_fasta), file(aln), file(ref_abund) from humann_faa_for_acc.join(humann2_ref_aln).join(ref_clustered_abund_humann2)
    val method_label from "HUMAnN2"

    output:
    file "${method_label}.${ix}.accuracy.tsv"

    script:
    template "calculate_gene_accuracy.sh"
}