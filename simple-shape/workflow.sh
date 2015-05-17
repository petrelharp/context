#!/bin/bash

set -eu
set -o pipefail

BASEDIR="testing"
GMDIR="genmatrices"
MODEL="shape-model-random-values.json"

mkdir -p $BASEDIR

echo "simulate up some sequence for testing"
SIMFILE="$BASEDIR/sim.RData"
SIMGENMAT="$GMDIR/sim-genmatrix-$(echo $MODEL | sed -e 's/.json//').RData"  # this takes a WHILE, so let's save it for future use
Rscript ../sim-seq.R -c $MODEL -t .3 -s 1000000 -o $SIMFILE -m $SIMGENMAT

echo "and count the Tmers"
LONGWIN=5
SHORTWIN=2
LEFTWIN=1
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
