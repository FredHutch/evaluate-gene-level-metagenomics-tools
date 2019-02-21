#!/bin/bash

set -e

nextflow \
    run \
    evaluate_gene_detection.nf \
    --genome_list test_data/test_genomes.csv.gz \
    --outdir nf_outputs \
    --num_genomes 2 \
    --mean_depth 50 \
    --max_depth 5 \
    --log_std 1 \
    -with-docker "ubuntu:16.04" \
    -resume
