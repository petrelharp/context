#!/bin/bash

set -eu
set -o pipefail

echo "simulate up some sequence for testing"
Rscript ../scripts/sim-seq.R -c cpg-model-constrained.json -t .1 -s 1000000 -d testsim -o test-cpg-01.RData
Rscript ../scripts/sim-seq.R -c neg-cpg-model-constrained.json -t .1 -s 1000000 -d testsim -o test-cpg-02.RData

echo "and count the Tmers"
Rscript ../scripts/count-seq.R -i testsim/test-cpg-01.RData -w 4 -s 2 -l 1 -o testsim/test-cpg-01-counts-3.RData
Rscript ../scripts/count-seq.R -i testsim/test-cpg-02.RData -w 4 -s 2 -l 1 -o testsim/test-cpg-02-counts-3.RData

echo "precompute generator matrices:"
Rscript ../scripts/make-genmat.R -c cpg-model-constrained.json -w 4 -o genmatrices/genmatrix-4-cpg.RData

echo "check simulated model matches expected"
../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd testsim/test-cpg-01.RData genmatrices/genmatrix-4-cpg.RData
../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd testsim/test-cpg-02.RData genmatrices/genmatrix-4-cpg.RData

echo "fit a model"
Rscript ../scripts/fit-model.R -t .1 -i testsim/test-cpg-01-counts-3.RData -l 1 -m genmatrices/genmatrix-4-cpg.RData -o testsim/test-cpg-01-fit.RData -c cpg-model-constrained.json
Rscript ../scripts/fit-model.R -t .1 -i testsim/test-cpg-02-counts-3.RData -l 1 -m genmatrices/genmatrix-4-cpg.RData -o testsim/test-cpg-02-fit.RData -c cpg-model-constrained.json

echo "do mcmc"
Rscript ../scripts/mcmc-model.R -t .1 -i testsim/test-cpg-01-fit.RData -o testsim/test-cpg-01-mcmc-1.RData -c cpg-model-constrained.json --blen 3
Rscript ../scripts/mcmc-model.R -t .1 -i testsim/test-cpg-01-mcmc-1.RData -o testsim/test-cpg-01-mcmc-2.RData -c cpg-model-constrained.json --blen 100 -b 1000



echo "look at results: MLE"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd testsim/test-cpg-01-fit.RData testsim/test-cpg-01.RData
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd testsim/test-cpg-02-fit.RData testsim/test-cpg-02.RData
echo "look at results: MCMC"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd testsim/test-cpg-01-mcmc-2.RData testsim/test-cpg-01.RData 


# stop here for now
exit 0

echo "fit a model started from the truth"
Rscript ../scripts/fit-model.R -t .1 -i ${BASEDIR}$BASE-123456-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}.counts -l ${LEFTWIN} -m $GENMAT -j 54321 -c ${MODEL} -o ${BASEDIR}test-cpg-123456-test-start.RData

