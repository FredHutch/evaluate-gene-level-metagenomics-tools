#!/bin/bash

set -e

nextflow \
    -C "nextflow.local.config" \
    run \
    evaluate_gene_detection.nf \
    --genome_list test_data/test_genomes.csv.gz \
    --refdb genes.fastp.gz.rep.fasta \
    --outdir nf_outputs \
    --num_genomes 1 \
    --num_simulations 3 \
    --mean_depth 50 \
    --max_depth 10 \
    --log_std .1 \
    -with-docker "ubuntu:16.04" \
    -work-dir "work/" \
    -resume
