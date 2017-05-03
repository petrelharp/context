#!/bin/bash

set -eu
set -o pipefail

if [ $# -gt 0 ]
then
    NRUNS=$1
else
    echo "Usage:   $0 (number of runs)"
    exit 0
fi

# precompute generator matrices:
#   the `complete.json` file is the same as ising-model.json but without model fitting stuff
for k in 4 5 6
do
    if ! [ -f genmatrices/genmatrix-${k}-complete.RData ]
    then
        Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w $k
    fi
done

# this takes ~700 secs to simulate
SEQLEN=1000000
TLEN=1.0
# where to put stuff
BASEDIR=coverage-runs
# the model for inference
MODEL=ising-model.json
# number of MCMC batches of length 100 each
MCMCITER=10000

for N in $(seq $NRUNS)
do
    (
        DIR=${BASEDIR}/$(printf "%05g" $RANDOM);
        echo "Simulation $N, in $DIR,"
        mkdir -p $DIR
        # simulate up some sequence for testing;
        Rscript ../scripts/sim-seq.R -c $MODEL -t $TLEN -s $SEQLEN -d $DIR -o ising.RData;
        # and count the Tmers;
        Rscript ../scripts/count-seq.R -i $DIR/ising.RData -w 4 -s 2 -l 1;
        Rscript ../scripts/count-seq.R -i $DIR/ising.RData -w 5 -s 3 -l 1;
        Rscript ../scripts/count-seq.R -i $DIR/ising.RData -w 6 -s 2 -l 2;

        # fit the model;
        Rscript ../scripts/fit-model.R -c $MODEL -i $DIR/ising-4-root-2-tip-l1-shift0.counts -t $TLEN -m genmatrices/genmatrix-4-complete.RData -o $DIR/ising-fit-4-2-1.RData;
        Rscript ../scripts/fit-model.R -c $MODEL -i $DIR/ising-5-root-3-tip-l1-shift0.counts -t $TLEN -m genmatrices/genmatrix-5-complete.RData -o $DIR/ising-fit-5-3-1.RData;
        # Rscript ../scripts/fit-model.R -c $MODEL -i $DIR/ising-6-root-2-tip-l2-shift0.counts -t $TLEN -m genmatrices/genmatrix-6-complete.RData -o $DIR/ising-fit-6-2-2.RData;

        # and mcmc
        MCMCID=$RANDOM
        Rscript ../scripts/mcmc-model.R -i $DIR/ising-fit-4-2-1.RData -c ising-model.json -b $MCMCITER -l 1 -j $MCMCID
        # too long
        Rscript ../scripts/mcmc-model.R -i $DIR/ising-fit-5-3-1.RData -c ising-model.json -b $MCMCITER -j $MCMCID
        # waaay too long
        # Rscript ../scripts/mcmc-model.R -i $DIR/ising-fit-6-2-2.RData -c ising-model.json -b 1000 -j $MCMCID

        # gather results into json
        for RDATA in $DIR/ising-fit*.RData
        do
            ../scripts/gather-results.R --json -f $RDATA -s $DIR/ising.RData > ${RDATA%RData}.json
        done

    ) &
done

../scripts/collect-params-results.R $BASEDIR/*/ising-fit*.json > coverage_results.tsv

wait;
