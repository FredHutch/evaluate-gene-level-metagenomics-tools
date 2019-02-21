#!/bin/bash

set -e;

# Make sure the input file has data
(( $([[ "${fasta_in}" == *.gz ]] && gunzip -c "${fasta_in}" | wc -l || cat "${fasta_in}" | wc -l) > 1 ))

mmseqs createdb ${fasta_in} DB;
mmseqs cluster DB clu tmp --min-seq-id ${identity} --max-seqs 1000000 -c ${overlap} --threads ${cpu};
mmseqs createtsv DB DB clu ${base}.clusters.tsv;
mmseqs result2repseq DB clu ${base}.rep;
mmseqs result2flat DB DB ${base}.rep ${base}.rep.fasta --use-fasta-header

# Make sure the output file has data
(( $([[ "${base}.rep.fasta" == *.gz ]] && gunzip -c "${base}.rep.fasta" | wc -l || cat "${base}.rep.fasta" | wc -l) > 1 ))
(( $([[ "${base}.clusters.tsv" == *.gz ]] && gunzip -c "${base}.clusters.tsv" | wc -l || cat "${base}.clusters.tsv" | wc -l) > 1 ))
