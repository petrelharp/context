#!/bin/bash

set -eu
set -o pipefail

if [ $# -lt 1 ]
then
    echo "Usage: setup-sim-genmats.sh (prefix to .json config file)"
    exit 1
fi

MODEL=$1

BASEDIR="minimal"
GMDIR="genmatrices"
mkdir -p $GMDIR

SIMFILE="$BASEDIR/sim_${MODEL}.RData"
SIMGENMAT="$GMDIR/sim-${MODEL}-genmatrix.RData"  # this takes a WHILE, so let's save it for future use
echo "Simulating from ${MODEL} , saving to ${SIMFILE} , writing simulation genmatrix to ${SIMGENMAT} ."
Rscript ../scripts/sim-seq.R -c $MODEL.json -t .0001 -s 20 -o $SIMFILE -m $SIMGENMAT
