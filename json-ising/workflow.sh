#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c ising-model.json -t .1 -s 10000 -o sim-ising-123456.RData

# and count the Tmers
Rscript ../count-seq.R -i sim-ising-123456.RData -w 3 -s 1 -l 1

# precompute generator matrices:
#   width-{3,4,5,6} with all selection
Rscript ../flex-genmat.R -c genmatrices/complete.json -w 3

# fit a model
Rscript ../fit-model.R -i sim-ising-123456.3.1.l1.counts -l 1 -m genmatrices/genmatrix-3-complete.RData -j 54321

# mcmc also
Rscript ../mcmc-model.R -i sim-ising-123456-genmatrix-6-complete-54321.RData -c ising-prior.json -b 3 -j 2222
Rscript ../mcmc-model.R -i sim-ising-123456-genmatrix-6-complete-54321-mcmc-2222.RData -b 3 -j 3333


# OK, now do this on more Tmer sizes:
Rscript ../count-seq.R -i sim-ising-123456.RData -w 4 -s 2 -l 1
Rscript ../flex-genmat.R -c genmatrices/complete.json -w 4
Rscript ../fit-model.R -i sim-ising-123456.4.2.l1.counts -l 1 -m genmatrices/genmatrix-4-complete.RData -j 54321

Rscript ../count-seq.R -i sim-ising-123456.RData -w 5 -s 3 -l 1
Rscript ../flex-genmat.R -c genmatrices/complete.json -w 5
Rscript ../fit-model.R -i sim-ising-123456.5.3.l1.counts -l 1 -m genmatrices/genmatrix-5-complete.RData -j 54321

Rscript ../count-seq.R -i sim-ising-123456.RData -w 6 -s 2 -l 2
Rscript ../flex-genmat.R -c genmatrices/complete.json -w 6
Rscript ../fit-model.R -i sim-ising-123456.6.2.l2.counts -l 2 -m genmatrices/genmatrix-6-complete.RData -j 54321
