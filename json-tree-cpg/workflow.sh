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
Rscript ../count-seq.R -i $BASEDIR/sim-tree-cpg.RData -c sp2 -w $LONGWIN -s $SHORTWIN -l $LEFTWIN

# precompute generator matrices:
#   width-3
GENMAT=genmatrices/genmatrix-${LONGWIN}-cpg.RData
Rscript ../make-genmat.R -c tree-cpg-model.json -w ${LONGWIN} -o ${GENMAT} -n sp1

# fit a model
Rscript ../fit-model.R -i $BASEDIR/sim-tree-cpg-3-sp1-1-sp2-l1.counts -c tree-cpg-model.json -w ${LONGWIN} -l ${LEFTWIN} -j 001
Rscript ../fit-model.R -i $BASEDIR/sim-tree-cpg-3-sp2-1-sp1-l1.counts -c tree-cpg-model.json -w ${LONGWIN} -l ${LEFTWIN} -j 001

# compute residuals
Rscript ../compute-resids.R -i $BASE-123456-genmatrix-${LONGWIN}-cpg-54321.RData -w 3 -s 1 -l 1 -m ${GENMAT}
