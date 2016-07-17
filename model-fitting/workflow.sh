#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c cpg-model.json -t .1 -s 1000 -d simseqs -o sim-cpg-123456.RData
# and count the Tmers, with long windows of width 4 = 1+2+1 and short windows of length 2.
Rscript ../count-seq.R -i sim-cpg-123456.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-4 with all one-base transitions
Rscript ../make-genmat.R -c genmatrices/singlebase.json -w 4 &
#   width-5 with all one-base transitions
Rscript ../make-genmat.R -c genmatrices/singlebase.json -w 5 &
#   width-2 with all one-base transitions
Rscript ../make-genmat.R -c genmatrices/singlebase.json -w 2 &
#   width-4 with all two-base transitions
Rscript ../make-genmat.R -c genmatrices/dualbases.json -w 4

# fit a simple model: width 4 with all one-base transitions on 1+2+1 Tmers
Rscript ../fit-model.R -i sim-cpg-123456.4.2.l1.counts -m genmatrices/genmatrix-4-singlebase.RData -o simple-fit.RData
# and a more complex one: width 4 with all two-base transitions
Rscript ../fit-model.R -i sim-cpg-123456.4.2.l1.counts -m genmatrices/genmatrix-4-dualbases.RData -o complex-fit.RData

# first look at simple model residuals on 1+2+1 Tmers (what it was fit with)
Rscript ../compute-resids.R -i simple-fit.RData

# and on residuals on shorter Tmers: 2+1+1
Rscript ../compute-resids.R -i simple-fit.RData -w 4 -s 1 -l 2

# now, longer residuals (length 5)
#   first compute counts of 1+3+1 Tmers
Rscript ../count-seq.R -i sim-cpg-123456.RData -w 5 -s 3 -l 1
#   then get the residuals on 2+1+2 Tmers
Rscript ../compute-resids.R -i simple-fit.RData -m genmatrices/genmatrix-5-singlebase.RData -c sim-cpg-123456.5.3.l1.counts -w 5 -s 1 -l 2

# and do some MCMC (er, not much...)
Rscript ../mcmc-model.R -i simple-fit.RData -c simple-priors.json -b 2 -l 2 -j 13579
# continue this run
Rscript ../mcmc-model.R -i simple-fit-mcmc-13579.RData -c simple-priors.json -b 2 -l 2
