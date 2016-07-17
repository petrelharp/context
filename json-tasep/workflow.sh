#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
echo "simulating"
Rscript ../sim-seq.R -c tasep-model.json -t 0.1 -s 10000 -d simseqs -o sim-tasep-123456.RData

# and count the Tmers
echo "counting"
Rscript ../count-seq.R -i simseqs/sim-tasep-123456.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-4 with all selection
echo "making genmatrix, width 4"
Rscript ../make-genmat.R -c genmatrices/complete.json -w 4 --meanboundary 1

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd simseqs/sim-tasep-123456.RData genmatrices/genmatrix-4-complete.RData

echo "fitting, width 4"
Rscript ../fit-model.R -i simseqs/sim-tasep-123456-4-root-2-tip-l1-shift0.counts -m genmatrices/genmatrix-4-complete.RData -c tasep-model.json -t 0.1 -j 54321

echo "computing residuals"
Rscript ../compute-resids.R -i simseqs/sim-tasep-123456-4-root-2-tip-l1-shift0-genmatrix-4-complete-54321.RData

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd simseqs/sim-tasep-123456-4-root-2-tip-l1-shift0-genmatrix-4-complete-54321.RData simseqs/sim-tasep-123456.RData
