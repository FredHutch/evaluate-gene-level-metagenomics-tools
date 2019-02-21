#!/bin/bash

set -e

python << END

from collections import defaultdict
import gzip

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

# Read in the cluster membership
clusters = dict([
  line.rstrip("\n").split("\t")[:2][::-1]
  for line in gzip_safe_open("${groups}")
])

# Add up the abundances for each cluster
clust_abund = defaultdict(float)
for line in gzip_safe_open("${abund_in}"):
    if ',' not in line:
        continue
    member, abund = line.rstrip("\n").split(",", 2)[:2]
    assert member in clusters, member + " not in clusters"
    clust_abund[clusters[member]] += float(abund)

with open("${base}.clust.${identity}.csv", "wt") as fo:
    for k, v in clust_abund.items():
        fo.write(k + ',' + str(v) + '\n')

END