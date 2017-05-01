#!/bin/bash

set -eu
set -o pipefail

if [[ -e '/home/rcf-40/pralph/panfs/context/cpg' ]]
then
    cd /home/rcf-40/pralph/panfs/context/cpg
    source /home/rcf-40/pralph/cmb/bin/R-setup-usc.sh
fi

MODEL="cpg-plus-epsilon.json"
SUBDIR="sim_no_CpG"
TLEN=0.1

MODELNAME="${MODEL%.json}"
# precompute generator matrices:
mkdir -p genmatrices
for LONGWIN in 3 4 5 6
do
    GENMAT="genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.RData"
    if [[ ! -e $GENMAT ]]
    then
        Rscript ../scripts/make-genmat.R -c $MODEL -w ${LONGWIN} -o "genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.RData"
    fi
done

SEQLEN=1000000
# where to put stuff
BASEDIR=simseqs/model_selection


for N in $(seq 16)
do
    (
        DIR=$BASEDIR/sim-$(printf "%05g" $RANDOM);
        echo "Simulation $N, in $DIR";
        mkdir -p "$DIR";
        # simulate up some sequence for testing;
        Rscript ../scripts/sim-seq.R -c $MODEL -t $TLEN -s $SEQLEN -d "$DIR" -o "sim.RData";
        for LONGWIN in 3 4 5 6;
        do
            SHORTWIN=$(( LONGWIN/2 ));
            SHORTWIN=$(( SHORTWIN>2?2:SHORTWIN ));
            LEFTWIN=$(( (LONGWIN-SHORTWIN)/2 ));
            # and count the Tmers;
            Rscript ../scripts/count-seq.R -i $DIR/sim.RData -w $LONGWIN -s $SHORTWIN -l $LEFTWIN;
            # fit the model;
            FITFILE="$DIR/test-cpg-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.RData" -o $FITFILE;
            Rscript ../scripts/gather-results.R --fit $FITFILE --sim $DIR/sim.RData --outfile $DIR/fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.json --json 2>/dev/null ;
        done;
    ) &
done

wait;
