#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c cpg-model.json -t .1 -s 1000 -d simseqs -o sim-cpg-123456.RData
# and count the paired tuples, width 4 = 1+2+1
Rscript ../count-seq.R -i sim-cpg-123456.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-4 with all one-base transitions
Rscript ../flex-genmat.R -c genmatrices/singlebase.json -w 4
#   width-5 with all one-base transitions
Rscript ../flex-genmat.R -c genmatrices/singlebase.json -w 5
#   width-4 with all two-base transitions
Rscript ../flex-genmat.R -c genmatrices/dualbases.json -w 4

# fit a simple model
Rscript ../fit-model.R -i sim-cpg-123456.4.2.l1.counts -l 1 -m genmatrices/genmatrix-4-singlebase.RData -o simple-fit.RData
# and a more complex one
Rscript ../fit-model.R -i sim-cpg-123456.4.2.l1.counts -l 1 -m genmatrices/genmatrix-4-dualbases.RData -o complex-fit.RData

# first look at residuals, width 4 (what it was fit with)
Rscript ../compute-resids.R -i simple-fit.RData

# and, shorter residuals:
Rscript ../compute-resids.R -i simple-fit.RData -w 4 -s 1 -l 2 

# now, longer residuals (length 5)
#   first compute counts
Rscript ../count-seq.R -i sim-cpg-123456.RData -w 5 -s 3 -l 1
#   then get the residuals
Rscript ../compute-resids.R -i simple-fit.RData -m genmatrices/genmatrix-5-singlebase.RData -c sim-cpg-123456.5.3.l1.counts -w 5 -s 1 -l 2 

