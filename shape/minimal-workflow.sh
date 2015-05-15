#!/bin/bash

set -eu
set -o pipefail

BASEDIR="minimal"
MODEL="shape-model-random-values.json"

mkdir -p $BASEDIR

echo "simulate up some sequence for testing"
SIMFILE="$BASEDIR/sim.RData"
Rscript ../sim-seq.R -c $MODEL -t .01 -s 1000 -o $SIMFILE

echo "and count the Tmers"
LONGWIN=5
SHORTWIN=1
LEFTWIN=2
GENMAT="$BASEDIR/genmatrix-${LONGWIN}.RData"

COUNTFILE=$BASEDIR/sim-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.counts
ls $SIMFILE && \
    Rscript ../count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE

echo "precompute generator matrices:"
Rscript ../make-genmat.R -c $MODEL -w ${LONGWIN} -o ${GENMAT}

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd ${BASEDIR}/sim.RData ${GENMAT}

echo "fit a model"
FITFILE=$BASEDIR/fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData
ls ../fit-model.R $COUNTFILE $MODEL $GENMAT && \
    Rscript ../fit-model.R -i $COUNTFILE -t .01 --maxit 5 -c $MODEL -m $GENMAT -o $FITFILE
