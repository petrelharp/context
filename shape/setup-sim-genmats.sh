#!/bin/bash

set -eu
set -o pipefail

BASEDIR="minimal"
GMDIR="genmatrices"
mkdir -p $GMDIR

for MODEL in shape-model-MGW shape-model-all-variables 
do
    echo "Simulating from ${MODEL} ."
    SIMFILE="$BASEDIR/sim_${MODEL}.RData"
    SIMGENMAT="$GMDIR/sim-${MODEL}-genmatrix.RData"  # this takes a WHILE, so let's save it for future use
    Rscript ../sim-seq.R -c $MODEL -t .0001 -s 20 -o $SIMFILE -m $SIMGENMAT
done
