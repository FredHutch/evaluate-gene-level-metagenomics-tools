# Simulate Metagenomes

In order to evaluate the accuracy of gene detection methods, it is
necessary to have access to a set of metagenomes with known genes
and their abundances. 

This is a very naive method for simulating metagenomes that uses an
arbitrary number of genomes selected randomly from RefSeq, with
abundances simulated randomly from a log-normal distribution.

The actual read simulations will be done with ART, and the genes
present in each genome will be ascertained directly from the RefSeq
annotations.


### Input

  * Text file with contents of RefSeq (from [here]
  (https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/) e.g. )
  * Number of genomes to simulate
  * Mean depth of sequencing


### Output

  * FASTQ file of simulated reads
  * FASTP file of the protein-coding gene sequences
  * CSV file with the depth of sequencing for each gene


### References

O'Leary NA, Wright MW, Brister JR, Ciufo S, Haddad D, McVeigh R, Rajput B, Robbertse B, Smith-White B, Ako-Adjei D, Astashyn A, Badretdin A, Bao Y, Blinkova O, Brover V, Chetvernin V, Choi J, Cox E, Ermolaeva O, Farrell CM, Goldfarb T, Gupta T, Haft D, Hatcher E, Hlavina W, Joardar VS, Kodali VK, Li W, Maglott D, Masterson P, McGarvey KM, Murphy MR, O'Neill K, Pujar S, Rangwala SH, Rausch D, Riddick LD, Schoch C, Shkeda A, Storz SS, Sun H, Thibaud-Nissen F, Tolstoy I, Tully RE, Vatsan AR, Wallin C, Webb D, Wu W, Landrum MJ, Kimchi A, Tatusova T, DiCuccio M, Kitts P, Murphy TD, Pruitt KD. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. 2016 Jan 4;44(D1):D733-45

Tatusova T, DiCuccio M, Badretdin A, Chetvernin V, Nawrocki EP, Zaslavsky L, Lomsadze A, Pruitt KD, Borodovsky M, Ostell J. NCBI prokaryotic genome annotation pipeline. Nucleic Acids Res. 2016 Aug 19;44(14):6614-24

Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth. ART: a next-generation sequencing read simulator, Bioinformatics (2012) 28 (4): 593-594