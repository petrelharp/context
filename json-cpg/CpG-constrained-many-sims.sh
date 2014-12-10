#!/bin/bash

set -eu
set -o pipefail

cd /home/rcf-40/pralph/panfs/context/json-cpg
source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh

MODEL="cpg-model-CpG-constrained.json"
SUBDIR="sim_no_CpG"
TLEN=0.1

MODELNAME="$(echo $MODEL|sed -e 's/.json//').RData"
# precompute generator matrices:
mkdir -p genmatrices
for LONGWIN in 3 4 5 6
do
    Rscript ../make-genmat.R -c $MODEL -w ${LONGWIN} -o "genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.json"
done

for N in $(seq 16)
do
    (
        DIR=$(printf "%05g" $RANDOM);
        echo "Simulation $N, in $SUBDIR/$DIR";
        mkdir -p "$SUBDIR/sim-$DIR";
        # simulate up some sequence for testing;
        Rscript ../sim-seq.R -c $MODEL -t $TLEN -s 100000 -d "$SUBDIR/sim-$DIR" -o "sim.RData";
        for LONGWIN in 3 4 5 6
        do
            SHORTWIN=$(( LONGWIN/2 ))
            SHORTWIN=$(( SHORTWIN>2?2:SHORTWIN ))
            LEFTWIN=$(( (LONGWIN-SHORTWIN)/2 ))
            # and count the Tmers;
            Rscript ../count-seq.R -i $SUBDIR/sim-$DIR/sim.RData -w $LONGWIN -s $SHORTWIN -l $LEFTWIN;
            # fit the model;
            FITFILE="$SUBDIR/sim-$DIR/test-cpg-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../fit-model.R -c $MODEL -t $TLEN -i "$SUBDIR/sim-$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.RData" -o $FITFILE
            Rscript ../gather-results.R --fit $FITFILE --sim $SUBDIR/sim-${DIR}/sim.RData --outfile $SUBDIR/sim-$DIR/fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.json --json 2>/dev/null 
        done
    ) &
done

wait;

# after, run:
#   Rscript ../collect-many-sims.R $(find simseqs -name "*json") > many-sims-results.tsv

