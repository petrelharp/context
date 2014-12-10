#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
Rscript ../sim-seq.R -c ising-model.json -t .1 -s 100000 -d simseqs -o test-ising-123.RData

# and count the Tmers
Rscript ../count-seq.R -i simseqs/test-ising-123.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-3
Rscript ../make-genmat.R -c genmatrices/complete.json -w 4

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd simseqs/test-ising-123.RData genmatrices/genmatrix-4-complete.RData 

echo "fitting the model"
Rscript ../fit-model.R -i simseqs/test-ising-123-4-root-2-tip-l1-shift0.counts -m genmatrices/genmatrix-4-complete.RData -c ising-model.json -t .1 -j 456

echo "computing residuals"
Rscript ../compute-resids.R -i simseqs/test-ising-123-4-root-2-tip-l1-shift0-genmatrix-4-complete-456.RData

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd simseqs/test-ising-123-4-root-2-tip-l1-shift0-genmatrix-4-complete-456.RData simseqs/test-ising-123.RData

echo "mcmc also"
Rscript ../mcmc-model.R -i simseqs/test-ising-123-4-root-2-tip-l1-shift0-genmatrix-4-complete-456.RData -c ising-model.json -t .1 -b 1000 -j 789

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd simseqs/test-ising-123-4-root-2-tip-l1-shift0-genmatrix-4-complete-456-mcmc-789.RData simseqs/test-ising-123.RData
