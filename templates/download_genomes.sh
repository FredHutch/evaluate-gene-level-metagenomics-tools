#!/bin/bash

set -e;

ending="_genomic.fna.gz"
unzip_ending="_genomic.fna"

cat ${genome_abund_csv} | sed 's/,.*//' | while read url; do
    suffix="\$(echo \$url | sed 's/.*\\///')"

    # Fetch the genome
    wget \$url/\$suffix\$ending
    
    gunzip \$suffix\$ending

    # Annotate the genome
    prokka --prefix \$suffix --fast --quiet \$suffix\$unzip_ending
    mv \$suffix/* ./

    gzip \$suffix.fna
    gzip \$suffix.faa

done

tar cvf genomes.tar *.fna.gz *.faa.gz