#!/bin/bash

set -eu 
set -o pipefail

# simulate
BASEDIR=sim-01
Rscript ../sim-seq.R -c tree-cpg-model.json -s 10000 -d $BASEDIR -o sim-tree-cpg.RData

# and count the Tmers
LONGWIN=3
SHORTWIN=1
LEFTWIN=1

Rscript ../count-seq.R -i $BASEDIR/sim-tree-cpg.RData -c sp1 -w $LONGWIN -s $SHORTWIN -l $LEFTWIN

# precompute generator matrices:
#   width-3
GENMAT=genmatrix-${LONGWIN}-cpg.RData
Rscript ../make-genmat.R -c $MODEL -w ${LONGWIN} -o ${GENMAT}

# fit a model
Rscript ../fit-model.R -i $BASE-123456.${LONGWIN}.${SHORTWIN}.l${LEFTWIN}.counts -l ${LEFTWIN} -m $GENMAT -j 54321

# compute residuals
Rscript ../compute-resids.R -i $BASE-123456-genmatrix-${LONGWIN}-cpg-54321.RData -w 3 -s 1 -l 1 -m ${GENMAT}
