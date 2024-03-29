#!/bin/bash

set -eu 
set -o pipefail


BASEDIR=test-sims
LONGWIN=4
SHORTWIN=2
LEFTWIN=1

# simulate
mkdir -p $BASEDIR
Rscript ../scripts/sim-seq.R -c tree-cpg-model.json -s 10000 -d $BASEDIR -o sim-tree-cpg.RData

# and count the Tmers
Rscript ../scripts/count-seq.R -i $BASEDIR/sim-tree-cpg.RData -c sp1 -w $LONGWIN -s $SHORTWIN -l $LEFTWIN --shift $LONGWIN &
Rscript ../scripts/count-seq.R -i $BASEDIR/sim-tree-cpg.RData -c sp2 -w $LONGWIN -s $SHORTWIN -l $LEFTWIN --shift $LONGWIN &

wait;

# precompute generator matrices:
#   width-3
mkdir -p genmatrices
Rscript ../scripts/make-genmat.R -c tree-cpg-model.json -w ${LONGWIN}

# fit a model
for SHIFT in $(seq $LONGWIN)
do
    Rscript ../scripts/fit-tree-model.R -i $BASEDIR/sim-tree-cpg-${LONGWIN}-sp1-${SHORTWIN}-sp2-l${LEFTWIN}-shift${LONGWIN}.counts.$SHIFT -c tree-cpg-model.json -j 001 &
    Rscript ../scripts/fit-tree-model.R -i $BASEDIR/sim-tree-cpg-${LONGWIN}-sp2-${SHORTWIN}-sp1-l${LEFTWIN}-shift${LONGWIN}.counts.$SHIFT -c tree-cpg-model.json -j 001 &
done

wait
