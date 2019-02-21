#!/bin/bash

set -e;

python << END
import pandas as pd

plass = ['${sep="', '" plass}']
metaspades = ['${sep="', '" metaspades}']
megahit = ['${sep="', '" megahit}']
famli = ['${sep="', '" famli}']
all_diamond = ['${sep="', '" all_diamond}']
unique_diamond = ['${sep="', '" unique_diamond}']

assert len(famli) == len(plass)
assert len(all_diamond) == len(plass)
assert len(unique_diamond) == len(plass)

output = []

for ix, fp_arr in enumerate(zip(
  plass,
  metaspades,
  megahit,
  famli,
  all_diamond,
  unique_diamond
)):
    for fp in fp_arr:
        df = pd.read_table(fp)
        df["shard"] = ix
        output.append(df)

output = pd.concat(output)
output.to_csv("results.csv", index=False)
END