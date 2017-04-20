#!/bin/bash

set -eu
set -o pipefail

if [ $# -lt 1 ]
then
    echo "Usage:
    ./test-it.sh (name of .json model, without the .json) [more model names]
    "
else
    MODELS=$(echo "$@" | sed -e 's/[.]json\>//g')
fi

LONGWIN=7
SHORTWIN=3
LEFTWIN=2

for MODEL in $MODELS
do

    BASEDIR="test-it"
    GMDIR="genmatrices"
    INITMODEL=${MODEL} # change this to not start at true values
    SEED=$(printf "%06.0f" $RANDOM)

    mkdir -p $BASEDIR

    MODELFILE=${MODEL}.json
    INITMODELFILE=${INITMODEL}.json
    SIMFILE="$BASEDIR/${MODEL}-${SEED}-sim.RData"
    SIMGENMAT="$GMDIR/${MODEL}-sim-genmatrix.RData"  # this takes a WHILE, so let's save it for future use

    GENMAT="$GMDIR/${MODEL}-genmatrix-${LONGWIN}.RData"
    COUNTFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.counts/")
    FITFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.fit.RData/")
    RESIDFILE=$(echo $SIMFILE | sed -e "s/.RData/-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.resid.tsv/")

    mkdir -p $BASEDIR
    mkdir -p $GMDIR

    echo "Simulating from ${MODEL} ."
    Rscript ../scripts/sim-seq.R -c $MODELFILE -t 0.7 -s 1000000 -m $SIMGENMAT -z $SEED -o $SIMFILE

    echo "counting tuples, recording to ${COUNTFILE}"
    ls $SIMFILE && \
        Rscript ../scripts/count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE

    echo "precompute generator matrices, saving to ${GENMAT}"
    Rscript ../scripts/make-genmat.R -c $MODELFILE -w ${LONGWIN} -o ${GENMAT}

    echo "check simulated model matches expected"
    ../scripts/templated-Rmd.sh ../scripts/check-sim.Rmd $SIMFILE ${GENMAT}

    echo "fit a model, saving to ${FITFILE}"
    ls $COUNTFILE $INITMODELFILE $GENMAT && \
        Rscript ../scripts/fit-model.R -i $COUNTFILE -t 1.0 --maxit 500 -c $INITMODELFILE -m $GENMAT -o $FITFILE

    echo "computing residuals"
    echo "saving residuals to ${RESIDFILE}"
    ls $FITFILE && \
        Rscript ../scripts/compute-resids.R -i $FITFILE -o $RESIDFILE

    echo "look at results"
    ../scripts/templated-Rmd.sh ../scripts/simulation.Rmd $FITFILE $SIMFILE

done
