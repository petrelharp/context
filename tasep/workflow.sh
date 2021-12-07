#!/bin/bash

set -eu
set -o pipefail

# simulate up some sequence for testing
echo "simulating"
Rscript ../scripts/sim-seq.R -c tasep-model.json -t 0.1 -s 10000 -d simseqs -o sim-tasep-123456.RData

# and count the Tmers
echo "counting"
Rscript ../scripts/count-seq.R -i simseqs/sim-tasep-123456.RData -w 4 -s 2 -l 1

# precompute generator matrices:
#   width-4 with all selection
echo "making genmatrix, width 4"
Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w 4 --meanboundary 1

echo "check simulated model matches expected"
../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd simseqs/sim-tasep-123456.RData genmatrices/genmatrix-4-complete.RData

echo "fitting, width 4"
Rscript ../scripts/fit-model.R -i simseqs/sim-tasep-123456-4-root-2-tip-l1-shift0.counts -m genmatrices/genmatrix-4-complete.RData -c tasep-model.json -t 0.1 -j 54321

echo "computing residuals"
Rscript ../scripts/compute-resids.R -i simseqs/sim-tasep-123456-4-root-2-tip-l1-shift0-genmatrix-4-complete-nonnegative-54321.RData

echo "look at results"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd simseqs/sim-tasep-123456-4-root-2-tip-l1-shift0-genmatrix-4-complete-nonnegative-54321.RData simseqs/sim-tasep-123456.RData
