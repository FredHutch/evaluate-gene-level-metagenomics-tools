#!/bin/bash

set -e;

ending="_genomic.fna.gz"
unzip_ending="_genomic.fna"
protein_ending="_protein.faa.gz"
protein_unzip_ending="_protein.faa"

cat ${genome_abund_csv} | sed 's/,.*//' | while read url; do
    suffix="\$(echo \$url | sed 's/.*\\///')"

    # Fetch the genome
    wget \$url/\$suffix\$ending
    
    gunzip \$suffix\$ending

    # Annotate the genome with Prokka
    prokka --prefix \$suffix --fast --quiet \$suffix\$unzip_ending
    mv \$suffix/* ./

    gzip \$suffix.fna
    gzip \$suffix.faa

done

tar cvf genomes.${random_seed}.tar *.fna.gz *.faa.gz