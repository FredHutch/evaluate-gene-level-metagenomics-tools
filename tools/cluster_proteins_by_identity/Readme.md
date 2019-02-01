# Cluster Proteins by Identity

### Approach

It is often useful to group together proteins that are highly similar. 
In the most obvious example, if a simulated dataset has two genomes with
the identical protein in them, it would be impossible to tell which of
the two is being detected. This utility takes care of that problem.

The first use case is to cluster a set of proteins by sequence identity.
The second use case is to cluster a set of proteins and then also combine
the abundance CSV so that the abundances are summed for clusters of
proteins.
