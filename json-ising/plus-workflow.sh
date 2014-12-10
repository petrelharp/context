#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c ising-model-plus.json -t .1 -s 100000 -o simplus/plus-ising-01.RData

# and count the Tmers
Rscript ../count-seq.R -i simplus/plus-ising-01.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-3
Rscript ../make-genmat.R -c ising-model-plus.json -w 4 -o genmatrices/genmatrix-4-plus.RData

echo "fitting the model"
Rscript ../fit-model.R -i simplus/plus-ising-01-4-root-2-tip-l1-shift0.counts -m genmatrices/genmatrix-4-plus.RData -j 456

echo "computing residuals"
Rscript ../compute-resids.R -i simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456.RData

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456.RData simplus/plus-ising-01.RData

echo "mcmc also"
Rscript ../mcmc-model.R -i simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456.RData -c ising-model-plus.json -b 1000 -j 789

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456-mcmc-789.RData simplus/plus-ising-01.RData
