#!/bin/bash

set -e

nextflow \
    run \
    evaluate_gene_detection.nf \
    --genome_list test_data/test_genomes.csv.gz \
    --outdir nf_outputs \
    -with-docker "ubuntu:16.04" \
    -resume
