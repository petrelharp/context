#!/bin/bash

set -eu
set -o pipefail

BASEDIR="testing"
GMDIR="genmatrices"
MODEL="shape-model-random-values"

mkdir -p $BASEDIR

echo "simulate up some sequence for testing"
SIMGENMAT="$GMDIR/sim-genmatrix-${MODEL}.RData"  # should be precomputed
SIMFILE="test-shape-01.RData"
Rscript ../sim-seq.R -c ${MODEL}.json -t .1 -s 1000000 -d $BASEDIR -o $SIMFILE -m $SIMGENMAT

echo "and count the Tmers"
LONGWIN=5
SHORTWIN=1
LEFTWIN=2
GENMAT="$BASEDIR/genmatrix-${LONGWIN}.RData"

echo "and count the Tmers"
COUNTFILE=$BASEDIR/sim-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.counts
ls $BASEDIR/$SIMFILE && \
    Rscript ../count-seq.R -i $BASEDIR/$SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE || (echo "count-seq failed"; exit 1)

echo "precompute generator matrices:"
Rscript ../make-genmat.R -c ${MODEL}.json -w ${LONGWIN} -o ${GENMAT} || (echo "make-genmat failed"; exit 1)

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd ${BASEDIR}/${SIMFILE} ${GENMAT} 

echo "fit a model"
FITFILE=$BASEDIR/fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData
ls ../fit-model.R $COUNTFILE $MODEL.json $GENMAT && \
    Rscript ../fit-model.R -i $COUNTFILE -t .01 --maxit 5 -c $MODEL.json -m $GENMAT -o $FITFILE || (echo "fit-model failed"; exit 1)

echo "look at results: MLE"
../templated-Rmd.sh ../simulation.Rmd $FITFILE $BASEDIR/$SIMFILE

echo "do mcmc"
MCMCFILE=$BASEDIR/mcmc-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}_01.RData
Rscript ../mcmc-model.R -t .1 -i $FITFILE -o $MCMCFILE -c ${MODEL}.json --blen 3
Rscript ../mcmc-model.R -t .1 -i $FITFILE -o $MCMCFILE -c ${MODEL}.json --blen 100 -b 1000


echo "look at results: MCMC"
../templated-Rmd.sh ../simulation.Rmd $MCMCFILE $BASEDIR/$SIMFILE


exit 0
