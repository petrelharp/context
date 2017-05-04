#!/bin/bash

set -eu
set -o pipefail

if [ $# -lt 1 ]
then
    echo "Usage: quicktest.sh (prefix to .json config file)"
    exit 1
fi

MODEL=$(echo $1 | sed -e 's/.json$//')
MODELFILE=${MODEL}.json

SEED=00023

BASEDIR="quicktest"
GMDIR="genmatrices"
mkdir -p $BASEDIR
mkdir -p $GMDIR

echo "Simulating from ${MODEL} ."
SIMGENMAT="$GMDIR/sim-${MODEL}-genmatrix.RData"  # this should be already done
SIMFILE="$BASEDIR/simseq-${MODEL}-seed-${SEED}.RData"
Rscript ../scripts/sim-seq.R -c $MODELFILE -t 1.0 -s 100 -m $SIMGENMAT -z $SEED -o $SIMFILE

LONGWIN=5
SHORTWIN=3
LEFTWIN=1

GENMAT="$GMDIR/${MODEL}-genmatrix-${LONGWIN}.RData"
COUNTFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.counts/")

echo "counting tuples, recording to ${COUNTFILE}"
ls $SIMFILE && \
    Rscript ../scripts/count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE

echo "precompute generator matrices, saving to ${GENMAT}"
if [ ! -f $GENMAT ]
then
    Rscript ../scripts/make-genmat.R -c $MODELFILE -w ${LONGWIN} -o ${GENMAT}
else
    echo " ... ${GENMAT} already exists, using that one."
fi


echo "check simulated model matches expected"
../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd $SIMFILE ${GENMAT}

echo "fit a model"
FITFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.fit.RData/")
ls $COUNTFILE $MODELFILE $GENMAT && \
    Rscript ../scripts/fit-model.R -i $COUNTFILE -t 1.0 --maxit 5 -c $MODELFILE -m $GENMAT -o $FITFILE

echo "computing residuals"
RESIDFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.resid.tsv/")
ls $FITFILE && \
    Rscript ../scripts/compute-resids.R -i $FITFILE -o $RESIDFILE

echo "look at results"
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd $FITFILE $SIMFILE
