#!/bin/sh

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c ising-model.json -t .1 -s 10000 -d simseqs -o sim-ising-123456.RData

# and count the paired tuples
Rscript ../count-seq.R -i sim-ising-123456.RData -w 3 -s 1 -l 1

# precompute generator matrices:
#   width-3 with all selection
Rscript ../flex-genmat.R -c genmatrices/complete.json -w 3

Rscript ../fit-model.R -i sim-ising-123456.counts -l 1 -m genmatrices/genmatrix-3-complete.RData -j 54321
