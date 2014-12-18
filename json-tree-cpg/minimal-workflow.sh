#!/bin/bash

set -eu 
set -o pipefail


BASEDIR=minimal
LONGWIN=2
SHORTWIN=1
LEFTWIN=1

# simulate
mkdir -p $BASEDIR
SIMFILE="$BASEDIR/sim.RData"
echo "Simulating to $SIMFILE ."
Rscript ../sim-seq.R -c tree-cpg-model.json -s 100 -o $SIMFILE

# and count the Tmers
COUNTSFILE="$BASEDIR/sim.counts"
echo "Counting to $COUNTSFILE ."
Rscript ../count-seq.R -i $SIMFILE -c sp1 -w $LONGWIN -s $SHORTWIN -l $LEFTWIN --shift $LONGWIN -o $COUNTSFILE 


# precompute generator matrices:
mkdir -p genmatrices
GENMAT="genmatrices/cpg-${LONGWIN}.RData"  # (specified in config file)
echo "Genmatrix to $GENMAT ."
Rscript ../make-genmat.R -c tree-cpg-model.json -w ${LONGWIN} 

# fit a model
SHIFT=1
FITFILE="$BASEDIR/fit.RData"
echo "Fitting to $FITFILE ."
Rscript ../fit-tree-model.R -i ${COUNTSFILE}.${SHIFT} -c tree-cpg-model.json -o $FITFILE --maxit 3

# save results as json
RESFILE="$BASEDIR/fit.json"
Rscript ../gather-results.R --fit $FITFILE --sim $SIMFILE --outfile $RESFILE --json
