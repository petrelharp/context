#!/bin/bash

set -eu
set -o pipefail

cd /home/rcf-40/pralph/panfs/context/tree-cpg
source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh

MODEL=tree-cpg-model.json
SEQLEN=100000
MAXIT=100
NJOBS=16

# precompute generator matrices:
mkdir -p genmatrices
for LONGWIN in 3 4 5 6
do
    GENMAT="genmatrices/cpg-${LONGWIN}.RData"  # (specified in config file)
    echo "Genmatrix to $GENMAT ."
    Rscript ../scripts/make-genmat.R -c ${MODEL} -w ${LONGWIN} 
done

for N in $(seq $NJOBS)
do
    (
        DIR=$(printf "%05g" $RANDOM);
        echo "Simulation $N, in simseqs/$DIR";
        mkdir -p simseqs/sim-$DIR
        # simulate up some sequence for testing;
        SIMFILE="simseqs/sim-$DIR/tree-cpg.RData"
        Rscript ../scripts/sim-seq.R -c tree-cpg-model.json -s $SEQLEN -o $SIMFILE
        # and count the Tmers;
        for LONGWIN in 3 4 5 6
        do
            SHORTWIN=$(( LONGWIN>5 ? LONGWIN-4 : LONGWIN-2 ))
            LEFTWIN=$(( (LONGWIN-SHORTWIN)/2 ))
            echo "Doing $LONGWIN, $SHORTWIN, $LEFTWIN now."
            COUNTSFILE="simseqs/sim-$DIR/tree-cpg-${LONGWIN}-${SHORTWIN}.counts"
            Rscript ../scripts/count-seq.R -i $SIMFILE -c sp1 -w $LONGWIN -s $SHORTWIN -l $LEFTWIN --shift $LONGWIN -o $COUNTSFILE;
            # fit the model;
            FITFILE="simseqs/sim-$DIR/fit-${LONGWIN}-${SHORTWIN}.RData"
            Rscript ../scripts/fit-tree-model.R -i ${COUNTSFILE}.1 -c tree-cpg-model.json -o $FITFILE --maxit $MAXIT
            RESFILE="simseqs/sim-$DIR/fit-${LONGWIN}-${SHORTWIN}.json"
            Rscript ../scripts/gather-results.R --fit $FITFILE --sim $SIMFILE --outfile $RESFILE --json 2>/dev/null 
        done
    ) &
done

wait;

# after, run:
#   Rscript ../scripts/collect-params-results.R $(find simseqs -name "*json") > many-sims-results.tsv
# and then
#   x <- read.table("many-sims-results.tsv",header=TRUE,check.names=FALSE)
#   layout(matrix(1:16,nrow=4))
#   for (k in 1:15) { hist(x[,15+k], xlim=range(x[,15+k],x[,k]), main=names(x)[15+k]); abline(v=unique(x[,k]),col='red'); hist(x[x$longwin==6,15+k],col=adjustcolor("blue",.5),add=TRUE) }

