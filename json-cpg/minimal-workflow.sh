#!/bin/bash

set -eu
set -o pipefail

BASE=test-cpg
BASEDIR=simseqs/
MODEL=cpg-model.json

echo "simulate up some sequence for testing"
Rscript ../sim-seq.R -c $MODEL -t .1 -s 1000 -d $BASEDIR -o $BASE-123456.RData

echo "and count the Tmers"
LONGWIN=3
SHORTWIN=1
LEFTWIN=1
GENMAT=genmatrices/genmatrix-${LONGWIN}-cpg.RData

Rscript ../count-seq.R -i ${BASEDIR}$BASE-123456.RData -w $LONGWIN -s $SHORTWIN -l $LEFTWIN

echo "precompute generator matrices:"
Rscript ../make-genmat.R -c $MODEL -w ${LONGWIN} -o ${GENMAT}

#echo "check simulated model matches expected"
#../templated-Rmd.sh ../testing-code/check-sim.Rmd ${BASEDIR}$BASE-123456.RData  ${GENMAT}

echo "fit a model"
ls ../fit-model.R ${BASEDIR}$BASE-123456-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}.counts $MODEL $GENMAT && \
Rscript ../fit-model.R -i ${BASEDIR}$BASE-123456-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}.counts -t .1 -c $MODEL -m $GENMAT -j 54321
