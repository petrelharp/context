#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c cpg-model.json -t .1 -s 10000 -d simseqs -o sim-cpg-123456.RData
# and count the paired tuples
Rscript ../count-seq.R -i sim-cpg-123456.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-5 with all one-base transitions
Rscript ../flex-genmat.R -c genmatrices/singlebase.json -w 4
#   width-5 with all two-base transitions
Rscript ../flex-genmat.R -c genmatrices/dualbases.json -w 4

# fit a simple model
Rscript ../fit-model.R -i sim-cpg-123456.counts -l 1 -m genmatrices/genmatrix-4-singlebase.RData -j 54321
# and a more complex one
Rscript ../fit-model.R -i sim-cpg-123456.counts -l 1 -m genmatrices/genmatrix-4-dualbases.RData -j 54321
