#!/bin/bash

set -eu
set -o pipefail

USAGE="Usage:\n run-analysis.sh (config file) (subdirectory) (long window) (short window) (length of time) (number of bases)\n"

if [ $# -lt 6 ]
then
    echo $USAGE
fi

MODEL=$1
BASEDIR=$2
LONGWIN=$3
SHORTWIN=$4
TLEN=$5
NSITES=$6

SUBDIR=$(printf "%06d" $RANDOM)
BASEDIR=$BASEDIR/sim-$SUBDIR

BASE=$(basename $MODEL | sed -e 's/.json$//')
LEFTWIN=$(( ($LONGWIN-$SHORTWIN)/2 ))

echo "Simulating some sequence from ${MODEL} ."
SIMFILE="${BASEDIR}/${BASE}-sim.RData"
Rscript ../scripts/sim-seq.R -c $MODEL -t $TLEN -s $NSITES -o $SIMFILE

GENMAT="genmatrices/genmatrix-${LONGWIN}-${BASE}.RData"
if [[ ! -e $GENMAT ]]
then
    echo "Making generator matrix $GENMAT ."
    mkdir -p genmatrices
    Rscript ../scripts/make-genmat.R -c $MODEL -w ${LONGWIN} -o $GENMAT
fi

COUNTFILE="$BASEDIR/$BASE-counts-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}-shift0.counts"
echo "Counting ${LONGWIN}/${SHORTWIN}/${LEFTWIN} Tmers, to $COUNTFILE ."
Rscript ../scripts/count-seq.R -i $SIMFILE -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTFILE

FITFILE="$BASEDIR/$BASE-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}-shift0.RData"
echo "Fitting a model, to $FITFILE ."
Rscript ../scripts/fit-model.R -i $COUNTFILE -m $GENMAT -c $MODEL -t $TLEN -o $FITFILE

echo "Look at results."
../scripts/templated-Rmd.sh ../scripts/simulation.Rmd $FITFILE $SIMFILE

## STOP HERE FOR NOW
exit 0

echo "mcmc also"
Rscript ../scripts/mcmc-model.R -i $FITFILE -c $MODEL -t $TLEN -b 1000 -j 3333
