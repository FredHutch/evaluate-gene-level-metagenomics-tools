#!/bin/bash

set -e;

# Make sure the input file has data
(( \$([[ "${fasta_in}" == *.gz ]] && gunzip -c "${fasta_in}" | wc -l || cat "${fasta_in}" | wc -l) > 1 ))

mmseqs createdb ${fasta_in} DB;
mmseqs cluster DB clu tmp --min-seq-id ${identity} --max-seqs 1000000 -c ${overlap} --threads 1;
mmseqs createtsv DB DB clu ${fasta_in}.clusters.tsv;
mmseqs result2repseq DB clu ${fasta_in}.rep;
mmseqs result2flat DB DB ${fasta_in}.rep ${fasta_in}.rep.fasta --use-fasta-header

# Make sure the output file has data
(( \$([[ "${fasta_in}.rep.fasta" == *.gz ]] && gunzip -c "${fasta_in}.rep.fasta" | wc -l || cat "${fasta_in}.rep.fasta" | wc -l) > 1 ))
(( \$([[ "${fasta_in}.clusters.tsv" == *.gz ]] && gunzip -c "${fasta_in}.clusters.tsv" | wc -l || cat "${fasta_in}.clusters.tsv" | wc -l) > 1 ))
