#!/bin/bash

set -eu
set -o pipefail

source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh

# precompute generator matrices:
Rscript ../make-genmat.R -c genmatrices/complete.json -w 3
Rscript ../make-genmat.R -c genmatrices/complete.json -w 4
Rscript ../make-genmat.R -c genmatrices/complete.json -w 5
Rscript ../make-genmat.R -c genmatrices/complete.json -w 6

for N in $(seq 16)
do
    (
        DIR=$RANDOM;
        echo "Simulation $N, in simseqs/$DIR";
        # simulate up some sequence for testing;
        Rscript ../sim-seq.R -c ising-model.json -t .1 -s 100000 -d simseqs/sim-$DIR -o test-ising.RData;
        # and count the Tmers;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 3 -s 1 -l 1;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 4 -s 2 -l 1;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 5 -s 3 -l 1;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-ising.RData -w 6 -s 2 -l 2;
        # fit the model;
        Rscript ../fit-model.R -i simseqs/sim-$DIR/test-ising-3-root-1-tip-l1-shift0.counts -m genmatrices/genmatrix-3-complete.RData -o simseqs/sim-$DIR/test-ising-fit-3-1-1.RData;
        Rscript ../fit-model.R -i simseqs/sim-$DIR/test-ising-4-root-2-tip-l1-shift0.counts -m genmatrices/genmatrix-4-complete.RData -o simseqs/sim-$DIR/test-ising-fit-4-2-1.RData;
        Rscript ../fit-model.R -i simseqs/sim-$DIR/test-ising-5-root-3-tip-l1-shift0.counts -m genmatrices/genmatrix-5-complete.RData -o simseqs/sim-$DIR/test-ising-fit-5-3-1.RData;
        Rscript ../fit-model.R -i simseqs/sim-$DIR/test-ising-6-root-2-tip-l2-shift0.counts -m genmatrices/genmatrix-6-complete.RData -o simseqs/sim-$DIR/test-ising-fit-6-2-2.RData;
        # look at results;
        ../templated-Rmd.sh ../simulation.Rmd simseqs/sim-$DIR/test-ising-fit-3-1-1.RData simseqs/sim-$DIR/test-ising.RData;
        ../templated-Rmd.sh ../simulation.Rmd simseqs/sim-$DIR/test-ising-fit-4-2-1.RData simseqs/sim-$DIR/test-ising.RData;
        ../templated-Rmd.sh ../simulation.Rmd simseqs/sim-$DIR/test-ising-fit-5-3-1.RData simseqs/sim-$DIR/test-ising.RData;
        ../templated-Rmd.sh ../simulation.Rmd simseqs/sim-$DIR/test-ising-fit-6-2-2.RData simseqs/sim-$DIR/test-ising.RData;
    ) &
done

wait;

