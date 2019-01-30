# Calculate Gene Accuracy


### Approach

Each sample was simulated from a set of genomes, where each genome contains a set of
genes and each genome is simulated with a given sequencing depth. By extension, each
sample has a set of genes that are actually present, and a sequencing depth for those
genes. A file summarizing the abundance of all of the genes from the simulation is
one of the inputs for this process.

For a given gene detection method, all of the detected genes have been aligned against
the set of genes that are actually present in the sample (from the simulation). That
file will be used to determine how accurately the method was able to detect the genes.

When analyzing the set of genes detected by a method, we ideally want to detect every
gene, with just a single gene reported for every gene that is present. In addition to
having the case where a gene is not detected (false positive), you may also have the
case where a method reports multiple genes as being present, when in truth there was
just one. We describe those as being "duplicate" detections. 

Operationally, "duplicate" detections are measured as the number of detected genes
that are not the top hit (most similar gene) for the truly present gene that is its
top hit. In other words, we take the protein-protein alignment table and sort by
alignment score, then we filter to only keep the top scoring alignment for every 
detected protein, and then we mark the top scoring alignment for every protein that
is truly present as the "true positive", with additional alignments as "duplicates."
In contrast, a "false positive" is a detected gene that does not align to the set
of genes that are truly present


### Metrics

  * Sensitivity: TP / (TP + FN)
  * Precision: TP / (TP + FP)
  * Uniqueness: TP / (TP + DUP)

Each metric can be calculated for the mock communities overall, or for the set of
genes that are simulated with a depth falling within a certain range (to get at the
question of limit of detection).

### Inputs

  * Alignment file comparing each detected gene to each gene that is actually present
  * TSV with the sequencing depth of each gene that is actually present in the simulation
  * String describing the detection method being evaluated

### Outputs

JSON file with the set of accuracy metrics, broken up by sequencing depth, e.g.:

| method | depth | sensitivity | precision | uniqueness |
| ------ | ----- | ----------- | --------- | ---------- |
| famli  | all   | 0.8         | 0.9       | 0.95       |
| famli  | 0.1   | 0.5         | 0.8       | 0.9        |
| famli  | 0.7   | 0.8         | 0.85      | 0.925      |
| famli  | 1.2   | 0.99        | 0.95      | 0.8        |
| plass  | all   | 0.9         | 0.8       | 0.99       |
| plass  | 0.1   | 0.6         | 0.7       | 0.8        |
| plass  | 0.7   | 0.9         | 0.75      | 0.905      |
| plass  | 1.2   | 0.95        | 0.85      | 0.8        |
