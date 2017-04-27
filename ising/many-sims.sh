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
for k in 3 4 5 6
do
    if ! [ -f genmatrices/genmatrix-${k}-complete.RData ]
    then
        Rscript ../scripts/make-genmat.R -c genmatrices/complete.json -w $k
    fi
done

# this takes ~200 secs to simulate
SEQLEN=1000000

for N in $(seq $NRUNS)
do
    (
        DIR=$(printf "%05g" $RANDOM);
        echo "Simulation $N, in simseqs/$DIR";
        mkdir -p simseqs/sim-$DIR
        # simulate up some sequence for testing;
        Rscript ../scripts/sim-seq.R -c ising-model.json -t .1 -s $SEQLEN -d simseqs/sim-$DIR -o test-ising.RData;
        # and count the Tmers;
        Rscript ../scripts/count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 3 -s 1 -l 1;
        Rscript ../scripts/count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 4 -s 2 -l 1;
        Rscript ../scripts/count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 5 -s 3 -l 1;
        Rscript ../scripts/count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 6 -s 2 -l 2;
        # fit the model;
        Rscript ../scripts/fit-model.R -i simseqs/sim-$DIR/test-ising-3-root-1-tip-l1-shift0.counts -t .1 -m genmatrices/genmatrix-3-complete.RData -o simseqs/sim-$DIR/test-ising-fit-3-1-1.RData;
        Rscript ../scripts/fit-model.R -i simseqs/sim-$DIR/test-ising-4-root-2-tip-l1-shift0.counts -t .1 -m genmatrices/genmatrix-4-complete.RData -o simseqs/sim-$DIR/test-ising-fit-4-2-1.RData;
        Rscript ../scripts/fit-model.R -i simseqs/sim-$DIR/test-ising-5-root-3-tip-l1-shift0.counts -t .1 -m genmatrices/genmatrix-5-complete.RData -o simseqs/sim-$DIR/test-ising-fit-5-3-1.RData;
        Rscript ../scripts/fit-model.R -i simseqs/sim-$DIR/test-ising-6-root-2-tip-l2-shift0.counts -t .1 -m genmatrices/genmatrix-6-complete.RData -o simseqs/sim-$DIR/test-ising-fit-6-2-2.RData;
    ) &
done

wait;

