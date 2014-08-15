#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing

Rscript ../sim-seq.R -c tasep-model.json -t 0.1 -s 10000 -d simseqs -o sim-tasep-123456.RData

# and count the Tmers
Rscript ../count-seq.R -i sim-tasep-123456.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-3 with all selection
Rscript ../flex-genmat.R -c genmatrices/complete.json -w 3

Rscript ../fit-model.R -i sim-tasep-123456.3.1.l1.counts -l 1 -m genmatrices/genmatrix-3-complete.RData -j 54321


