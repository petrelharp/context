#!/bin/bash

set -eu
set -o pipefail

echo "simulate up some sequence for testing"
Rscript ../sim-seq.R -c cpg-model-constrained.json -t .1 -s 1000000 -d testsim -o test-cpg-01.RData

echo "and count the Tmers"
Rscript ../count-seq.R -i testsim/test-cpg-01.RData -w 3 -s 1 -l 1 -o testsim/test-cpg-01-counts-3.RData

echo "precompute generator matrices:"
Rscript ../make-genmat.R -c cpg-model-constrained.json -w 3 -o genmatrices/genmatrix-3-cpg.RData

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd testsim/test-cpg-01.RData genmatrices/genmatrix-3-cpg.RData

echo "fit a model"
Rscript ../fit-model.R -i testsim/test-cpg-01-counts-3.RData -l 1 -m genmatrices/genmatrix-3-cpg.RData -o testsim/test-cpg-01-fit.RData -c cpg-model-constrained.json

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd testsim/test-cpg-01-fit.RData testsim/test-cpg-01.RData
