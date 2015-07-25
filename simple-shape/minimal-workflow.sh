#!/bin/bash

set -eu
set -o pipefail

BASEDIR="minimal"
GMDIR="genmatrices"
MODEL="simple-shape-model"
INITMODEL="simple-shape-model-init" # don't start at true values
SEED=$(printf "%06.0f" $RANDOM)

mkdir -p $BASEDIR

MODELFILE=${MODEL}.json
INITMODELFILE=${INITMODEL}.json
SIMFILE="$BASEDIR/sim.RData"
SIMGENMAT="$GMDIR/sim-genmatrix.RData"  # this takes a WHILE, so let's save it for future use

LONGWIN=7
SHORTWIN=3
LEFTWIN=2
GENMAT="$BASEDIR/genmatrix-${LONGWIN}.RData"

GENMAT="$GMDIR/${MODEL}-genmatrix-${LONGWIN}.RData"
COUNTFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.counts/")
FITFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.fit.RData/")
RESIDFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.resid.tsv/")

mkdir -p $BASEDIR
mkdir -p $GMDIR

echo "Simulating from ${MODEL} ."
Rscript ../sim-seq.R -c $MODELFILE -t 0.7 -s 1000000 -m $SIMGENMAT -z $SEED -o $SIMFILE

echo "counting tuples, recording to ${COUNTFILE}"
ls $SIMFILE && \
    Rscript ../count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE

echo "precompute generator matrices, saving to ${GENMAT}"
if [ ! -f $GENMAT ]
then
    Rscript ../make-genmat.R -c $MODELFILE -w ${LONGWIN} -o ${GENMAT}
else
    echo " ... ${GENMAT} already exists, using that one."
fi

echo "check simulated model matches expected"
../templated-Rmd.sh ../testing-code/check-sim.Rmd $SIMFILE ${GENMAT}

echo "fit a model, saving to ${FITFILE}"
ls $COUNTFILE $INITMODELFILE $GENMAT && \
    Rscript ../fit-model.R -i $COUNTFILE -t 1.0 --maxit 500 -c $INITMODELFILE -m $GENMAT -o $FITFILE

echo "computing residuals"
echo "saving residuals to ${RESIDFILE}"
ls $FITFILE && \
    Rscript ../compute-resids.R -i $FITFILE -o $RESIDFILE

echo "look at results"
../templated-Rmd.sh ../simulation.Rmd $FITFILE $SIMFILE

