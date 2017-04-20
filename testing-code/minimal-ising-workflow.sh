#!/bin/bash

set -eu
set -o pipefail

BASEDIR="minimal"
MODEL="../ising/ising-model.json"

mkdir -p $BASEDIR

echo "simulate up some sequence for testing"
SIMFILE="$BASEDIR/sim.RData"
Rscript ../scripts/sim-seq.R -c $MODEL -t .01 -s 1000 -o $SIMFILE

echo "and count the Tmers"
LONGWIN=2
SHORTWIN=1
LEFTWIN=1
GENMAT="$BASEDIR/genmatrix-${LONGWIN}.RData"

COUNTFILE=$BASEDIR/sim-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.counts
ls $SIMFILE && \
    Rscript ../scripts/count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE

echo "precompute generator matrices:"
Rscript ../scripts/make-genmat.R -c $MODEL -w ${LONGWIN} -o ${GENMAT}

#echo "check simulated model matches expected"
#../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd ${BASEDIR}$BASE-123456.RData  ${GENMAT}

echo "fit a model"
FITFILE=$BASEDIR/fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData
ls ../scripts/fit-model.R $COUNTFILE $MODEL $GENMAT && \
    Rscript ../scripts/fit-model.R -i $COUNTFILE -t .01 --maxit 5 -c $MODEL -m $GENMAT -o $FITFILE

echo "and look at output"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd $FITFILE $SIMFILE


exit 0;

# THIS IS SLIGHTLY LESS MINIMAL:

echo "now do it with shifts: counting again"
Rscript ../scripts/count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN --shift 2 -o $COUNTFILE

echo "and fitting"
FITFILE=$BASEDIR/fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}-s1.RData
ls ${COUNTFILE}.1 &&
    Rscript ../scripts/fit-model.R -i ${COUNTFILE}.1 -t .01 --maxit 5 -c $MODEL -m $GENMAT -o $FITFILE


