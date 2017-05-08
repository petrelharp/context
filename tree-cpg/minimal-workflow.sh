#!/bin/bash

set -eu 
set -o pipefail


BASEDIR=minimal
LONGWIN=3
SHORTWIN=2
LEFTWIN=1

# simulate
mkdir -p $BASEDIR
SIMFILE="$BASEDIR/sim.RData"
echo "Simulating to $SIMFILE ."
Rscript ../scripts/sim-seq.R -c tree-cpg-model.json -s 100 -o $SIMFILE

# and count the Tmers
COUNTSFILE="$BASEDIR/sim.counts"
echo "Counting to $COUNTSFILE ."
Rscript ../scripts/count-seq.R -i $SIMFILE -c sp1 -w $LONGWIN -s $SHORTWIN -l $LEFTWIN -o $COUNTSFILE 


# precompute generator matrices:
mkdir -p genmatrices
GENMAT="genmatrices/cpg-${LONGWIN}.RData"  # (specified in config file)
echo "Genmatrix to $GENMAT ."
Rscript ../scripts/make-genmat.R -c tree-cpg-model.json -w ${LONGWIN} 

# fit a model
FITFILE="$BASEDIR/fit.RData"
echo "Fitting to $FITFILE ."
Rscript ../scripts/fit-tree-model.R -i ${COUNTSFILE} -c tree-cpg-model.json -o $FITFILE --maxit 3

# save results as json
RESFILE="$BASEDIR/fit.json"
Rscript ../scripts/gather-results.R --fit $FITFILE --sim $SIMFILE --outfile $RESFILE --json

# compute residuals
RESIDFILE=$"$BASEDIR/resids-2-1-1.tsv"
Rscript ../scripts/compute-resids.R -i $FITFILE -o $RESIDFILE -w 2 -s 1 -l 1 --pretty 

