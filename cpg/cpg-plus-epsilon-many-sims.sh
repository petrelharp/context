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

MODEL="cpg-plus-epsilon.json"

BASEMODEL="base-model.json"
MODEL2="cpg-plus-2mers.json"
MODEL3="cpg-plus-3mers.json"

SUBDIR="simseqs/cpg-plus"
TLEN=0.3

# precompute generator matrices:
mkdir -p genmatrices
for THISMODEL in $MODEL $BASEMODEL $MODEL2 $MODEL3
do
    MODELNAME="${THISMODEL%.json}"
    for LONGWIN in 4 5 6
    do
        GENMAT="genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.RData"
        if [[ ! -e $GENMAT ]]
        then
            Rscript ../scripts/make-genmat.R -c $THISMODEL -w ${LONGWIN} -o "genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.RData"
        fi
    done
done

SEQLEN=1000000
# where to put stuff
BASEDIR=simseqs/model_selection


# for N in $(seq $NRUNS)
# do
#     (
#        DIR=$BASEDIR/$(printf "%05g" $RANDOM);

DIR=$BASEDIR/12474
        echo "Simulation $N, in $DIR";
        mkdir -p "$DIR";
        # simulate up some sequence for testing;
        if [ ! -f $DIR/sim.RData ]
        then
            Rscript ../scripts/sim-seq.R -c $MODEL -t $TLEN -s $SEQLEN -d "$DIR" -o "sim.RData";
        fi

        for LONGWIN in 4 5;
        do
            # this gets 1-2-1 and 1-3-1
            SHORTWIN=$(( (LONGWIN+1)/2 ));
            LEFTWIN=$(( (LONGWIN-SHORTWIN)/2 ));

            # count the Tmers;
            Rscript ../scripts/count-seq.R -i $DIR/sim.RData -w $LONGWIN -s $SHORTWIN -l $LEFTWIN;

            # fit the model;
            MODELNAME=${BASEMODEL%.json}
            FITFILE="$DIR/base-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $BASEMODEL -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME}.RData" -o $FITFILE;
            Rscript ../scripts/gather-results.R --fit $FITFILE --sim $DIR/sim.RData --json > ${FITFILE%RData}json;

            # now with all 2mers
            MODELNAME2="${MODEL2%.json}"
            FITFILE="$DIR/cpg-2mer-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL2 -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME2}.RData" -o $FITFILE;
            Rscript ../scripts/gather-results.R --fit $FITFILE --sim $DIR/sim.RData --json;
        done;
#     ) &
# done

../scripts/collect-params-results.R $DIR/cpg-fit*.json > cpg-plus-epsilon_results.tsv

wait;
