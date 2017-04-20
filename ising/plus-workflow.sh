#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
mkdir -p simplus
Rscript ../scripts/sim-seq.R -c ising-model-plus.json -t .1 -s 100000 -o simplus/plus-ising-01.RData

# and count the Tmers
Rscript ../scripts/count-seq.R -i simplus/plus-ising-01.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-3
Rscript ../scripts/make-genmat.R -c ising-model-plus.json -w 4 -o genmatrices/genmatrix-4-plus.RData

echo "fitting the model"
Rscript ../scripts/fit-model.R -i simplus/plus-ising-01-4-root-2-tip-l1-shift0.counts -t .1 -m genmatrices/genmatrix-4-plus.RData -j 456 -c ising-model-plus.json

echo "look at results"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456.RData simplus/plus-ising-01.RData

echo "mcmc also"
Rscript ../scripts/mcmc-model.R -i simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456.RData -c ising-model-plus.json -t .1 -b 10
Rscript ../scripts/mcmc-model.R -i simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456.RData -c ising-model-plus.json -t .1 -b 1000 -j 789

echo "look at results"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simplus/plus-ising-01-4-root-2-tip-l1-shift0-genmatrix-4-plus-456-mcmc-789.RData simplus/plus-ising-01.RData
