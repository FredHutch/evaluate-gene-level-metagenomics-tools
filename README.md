# Reproducible Evaluation of Gene-Level Metagenomics Tools

Samuel S. Minot, PhD & Jonathan L. Golob, MD/PhD

### Purpose

To evaluate the performance of tools used for gene-level metagenomics, in other
words, tools that detect and quantify genes from complex mixtures of organisms.


### Background

Tools for detecting and quantifying genes from metagenomes (genomic data generated
from mixtures of organisms) typically fall into two categories:

  * Assemblers: Tools that reconstruct genomes or genes directly from WGS data
  * Mappers: Tools that compare WGS data to a reference database of genes


### Approach

Our overall approach is to:

  1) Simulate a set of metagenomes, in which it is known exactly what genes are present
  2) Run a set of tools over those simulated WGS datasets
  3) Evaluate the ability of each tool to efficiently detect each gene, as a function of sequencing depth


### Implementation

We hope to achieve reproducibility by building this evaluation with Cromwell, 
a workflow management system that coordinates the execution of commands within
Docker containers and the exchange of data objects between those tasks. More 
details on Cromwell [here](http://cromwell.readthedocs.io).


### Other Resources

Many of the tools used in this evaluation are written as Cromwell tasks within the ongoing
project for [reproducible workflows at Fred Hutch](https://github.com/FredHutch/reproducible-workflows).

Useful tool for visualizing Cromwell workflows: [http://pb.opensource.epam.com/](http://pb.opensource.epam.com/)
