# Check Unique Fasta Headers

When aligning sequences, it is important to have unique headers for each sequence
record in the FASTA file. This prevents the confusion of having two records align
to different things, but having those be indistinguishable in the alignment
report. Another problem that this tool addresses is when there are whitespaces
in the headers, which are interpreted differently by various tools.

This just makes sure that the input file has unique headers.

### Input

  * FASTA file

### Output

  * FASTA file
