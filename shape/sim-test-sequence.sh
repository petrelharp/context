#!/bin/bash

set -eu
set -o pipefail

if [ $# -lt 1 ]
then
    echo "Usage: sim-test-sequence.sh (prefix to .json config file)"
    exit 1
fi

MODEL=$(echo $1 | sed -e 's/.json$//')
MODELFILE=${MODEL}.json

SEED=$(printf "%06.0f" $RANDOM)

BASEDIR="testseq"
GMDIR="genmatrices"
mkdir -p $BASEDIR
mkdir -p $GMDIR

echo "Simulating from ${MODEL} ."
SIMGENMAT="$GMDIR/sim-${MODEL}-genmatrix.RData"  # this should be already done
SIMFILE="$BASEDIR/simseq-${MODEL}-seed-${SEED}.RData"
Rscript ../sim-seq.R -c $MODELFILE -t .2 -s 10000 -m $SIMGENMAT -z $SEED -o $SIMFILE

LONGWIN=5
SHORTWIN=1
LEFTWIN=2

GENMAT="$GMDIR/${MODEL}-genmatrix-${LONGWIN}.RData"
COUNTFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.counts/")

echo "counting tuples:"
ls $SIMFILE && \
    Rscript ../count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE

echo "precompute generator matrices:"
Rscript ../make-genmat.R -c $MODELFILE -w ${LONGWIN} -o ${GENMAT}

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd $SIMFILE ${GENMAT}

echo "fit a model"
FITFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.fit.RData/")
ls $COUNTFILE $MODELFILE $GENMAT && \
    Rscript ../fit-model.R -i $COUNTFILE -t .2 --maxit 500 -c $MODELFILE -m $GENMAT -o $FITFILE

echo "computing residuals"
RESIDFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.resid.tsv/")
ls $FITFILE && \
    Rscript ../compute-resids.R -i $FITFILE -o $RESIDFILE

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd $FITFILE $SIMFILE
