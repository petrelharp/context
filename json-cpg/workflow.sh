#!/bin/bash

set -eu
set -o pipefail

BASE=test-cpg
BASEDIR=simseqs/
MODEL=cpg-model.json

echo "simulate up some sequence for testing"
Rscript ../sim-seq.R -c $MODEL -t .1 -s 1000000 -d $BASEDIR -o $BASE-123456.RData

echo "and count the Tmers"
LONGWIN=4
SHORTWIN=2
LEFTWIN=1
GENMAT=genmatrices/genmatrix-${LONGWIN}-cpg.RData

Rscript ../count-seq.R -i ${BASEDIR}$BASE-123456.RData -w $LONGWIN -s $SHORTWIN -l $LEFTWIN

echo "precompute generator matrices:"
Rscript ../make-genmat.R -c $MODEL -w ${LONGWIN} -o ${GENMAT}

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd ${BASEDIR}$BASE-123456.RData  ${GENMAT}

echo "fit a model"
Rscript ../fit-model.R -i ${BASEDIR}$BASE-123456-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}.counts -l ${LEFTWIN} -m $GENMAT -j 54321

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd ${BASEDIR}$BASE-123456-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-genmatrix-${LONGWIN}-cpg-54321.RData ${BASEDIR}$BASE-123456.RData

echo "compute residuals"
Rscript ../compute-resids.R -i ${BASEDIR}$BASE-123456-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-genmatrix-${LONGWIN}-cpg-54321.RData -w ${LONGWIN} -s 2 -l 1 -m ${GENMAT}
Rscript ../compute-resids.R -i ${BASEDIR}$BASE-123456-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-genmatrix-${LONGWIN}-cpg-54321.RData -w 2 -s 2 -l 0 -m ${GENMAT}

## STOP HERE FOR NOW
exit 0

echo "mcmc also"
Rscript ../mcmc-model.R -i $BASE-123456-genmatrix-3-cpg-54321.RData -c ${MODEL} -b 3 -j 1111
Rscript ../mcmc-model.R -i $BASE-123456-genmatrix-3-cpg-54321-mcmc-1111.RData -b 100 -j 2222
Rscript ../mcmc-model.R -i $BASE-123456-genmatrix-3-cpg-54321-mcmc-2222.RData -b 1000 -j 3333


echo "OK, now do this on more Tmer sizes:"
Rscript ../count-seq.R -i $BASE-123456.RData -w 4 -s 2 -l 1
Rscript ../make-genmat.R -c ${MODEL} -w 4
Rscript ../fit-model.R -i $BASE-123456.4.2.l1.counts -l 1 -m genmatrix-4-cpg-model.RData -j 54321

Rscript ../count-seq.R -i $BASE-123456.RData -w 5 -s 3 -l 1
Rscript ../make-genmat.R -c ${MODEL} -w 5
Rscript ../fit-model.R -i $BASE-123456.5.3.l1.counts -l 1 -m genmatrix-5-cpg-model.RData -j 54321

Rscript ../count-seq.R -i $BASE-123456.RData -w 6 -s 2 -l 2
Rscript ../make-genmat.R -c ${MODEL} -w 6
Rscript ../fit-model.R -i $BASE-123456.6.2.l2.counts -l 2 -m genmatrix-6-cpg-model.RData -j 54321
