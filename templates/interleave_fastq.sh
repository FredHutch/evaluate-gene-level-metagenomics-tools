#!/bin/bash

set -e;

tar xvf "${all_reads_tar}"

for r1 in *1.fq.gz; do

    r2="\$(echo \$r1 | sed 's/1.fq.gz/2.fq.gz/')"

    if [[ -s \$r2 ]]; then

        python << END
import gzip

with gzip.open("\$r1", "rt") as f1, gzip.open("\$r2", "rt") as f2, open("reads.fastq", "wa") as fo:
    while True:
        line = f1.readline()
        if line.strip() == "":
            break
        fo.write(line)
        
        for i in range(3):
            fo.write(f1.readline())
        
        for i in range(4):
            fo.write(f2.readline())
END

    rm "\$r1" "\$r2"

    fi

done

gzip reads.fastq
