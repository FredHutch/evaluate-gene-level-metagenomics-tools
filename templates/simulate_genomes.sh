#!/bin/bash

set -e;

tar xvf "${genome_tar}"

cat "${genome_abund_csv}" | tr ',' '\t' | while read url depth; do
    acc="\$(echo \$url | sed 's/.*\\///')"
    [[ -s \$acc.fna.gz ]] || ( echo "Can't find \$acc.fna.gz" && break )

    gunzip \$acc.fna.gz

    art_illumina --fcov \$depth --in \$acc.fna --len ${read_length} --mflen ${mflen} --sdev ${sdev} --noALN --out \$acc. --paired --quiet
    gzip *[12].fq

done

tar cvf all_reads.tar *.fq.gz