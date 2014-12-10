#!/bin/bash

set -eu
set -o pipefail

cd /home/rcf-40/pralph/panfs/context/json-cpg
source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh

MODEL=cpg-model.json

# precompute generator matrices:
mkdir -p genmatrices
for LONGWIN in 3 4 5 6
do
    Rscript ../make-genmat.R -c $MODEL -w ${LONGWIN} -o genmatrices/genmatrix-${LONGWIN}-cpg.RData
done

for N in $(seq 16)
do
    (
        DIR=$(printf "%05g" $RANDOM);
        echo "Simulation $N, in simseqs/$DIR";
        mkdir -p simseqs/sim-$DIR
        # simulate up some sequence for testing;
        Rscript ../sim-seq.R -c $MODEL -t .1 -s 100000 -d simseqs/sim-$DIR -o test-cpg.RData;
        # and count the Tmers;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-cpg.RData -w 3 -s 1 -l 1;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-cpg.RData -w 4 -s 2 -l 1;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-cpg.RData -w 5 -s 3 -l 1;
        Rscript ../count-seq.R -i simseqs/sim-$DIR/test-cpg.RData -w 6 -s 2 -l 2;

        # fit the model;
        Rscript ../fit-model.R -c $MODEL -t .1 -i simseqs/sim-$DIR/test-cpg-3-root-1-tip-l1-shift0.counts -m genmatrices/genmatrix-3-cpg.RData -o simseqs/sim-$DIR/test-cpg-fit-3-1-1.RData;
        Rscript ../gather-results.R --fit simseqs/sim-$DIR/test-cpg-fit-3-1-1.RData --sim simseqs/sim-${DIR}/test-cpg.RData --outfile simseqs/sim-$DIR/test-cpg-fit-3-1-1.json --json 2>/dev/null 

        Rscript ../fit-model.R -c $MODEL -t .1 -i simseqs/sim-$DIR/test-cpg-4-root-2-tip-l1-shift0.counts -m genmatrices/genmatrix-4-cpg.RData -o simseqs/sim-$DIR/test-cpg-fit-4-2-1.RData;
        Rscript ../gather-results.R --fit simseqs/sim-$DIR/test-cpg-fit-4-2-1.RData --sim simseqs/sim-${DIR}/test-cpg.RData --outfile simseqs/sim-$DIR/test-cpg-fit-4-2-1.json --json 2>/dev/null 

        Rscript ../fit-model.R -c $MODEL -t .1 -i simseqs/sim-$DIR/test-cpg-5-root-3-tip-l1-shift0.counts -m genmatrices/genmatrix-5-cpg.RData -o simseqs/sim-$DIR/test-cpg-fit-5-3-1.RData;
        Rscript ../gather-results.R --fit simseqs/sim-$DIR/test-cpg-fit-5-3-1.RData --sim simseqs/sim-${DIR}/test-cpg.RData --outfile simseqs/sim-$DIR/test-cpg-fit-5-3-1.json --json 2>/dev/null 

        Rscript ../fit-model.R -c $MODEL -t .1 -i simseqs/sim-$DIR/test-cpg-6-root-2-tip-l2-shift0.counts -m genmatrices/genmatrix-6-cpg.RData -o simseqs/sim-$DIR/test-cpg-fit-6-2-2.RData;
        Rscript ../gather-results.R --fit simseqs/sim-$DIR/test-cpg-fit-6-2-2.RData --sim simseqs/sim-${DIR}/test-cpg.RData --outfile simseqs/sim-$DIR/test-cpg-fit-6-2-2.json --json 2>/dev/null 
    ) &
done

wait;

# after, run:
#   Rscript ../collect-many-sims.R $(find simseqs -name "*json") > many-sims-results.tsv

