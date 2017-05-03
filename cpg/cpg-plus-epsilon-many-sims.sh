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
MODEL2="base-plus-2mers.json"
MODEL3="base-plus-3mers.json"

MODEL_A="cpg-only.json"
MODEL_C="cpg-plus-models/cpg-plus-1.json"
MODEL_D="cpg-plus-models/cpg-plus-2.json"

SUBDIR="simseqs/cpg-plus"
TLEN=0.3

# precompute generator matrices:
mkdir -p genmatrices
for THISMODEL in $MODEL $BASEMODEL $MODEL2 $MODEL3 $MODEL_A $MODEL_C $MODEL_D
do
    MODELNAME=$(basename ${THISMODEL%.json})
    for LONGWIN in 4 5
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

            # compute residuals
            ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
            ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 

            # now with all 2mers
            MODELNAME2="${MODEL2%.json}"
            FITFILE2="$DIR/2mer-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL2 -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME2}.RData" -o $FITFILE2;
            Rscript ../scripts/gather-results.R --fit $FITFILE2 --sim $DIR/sim.RData --json > ${FITFILE2%RData}json;

            # CLEAR FROM RESIDUALS: now with CpG
            MODELNAME_A=${MODEL_A%.json}
            FITFILE_A="$DIR/${MODELNAME_A}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL_A -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME_A}.RData" -o $FITFILE_A;
            Rscript ../scripts/gather-results.R --fit $FITFILE_A --sim $DIR/sim.RData --json > ${FITFILE_A%RData}json;
            ../scripts/compute-resids.R -i $FITFILE_A -o ${FITFILE_A%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 

            # AGAIN CLEAR FROM RESIDUALS
            MODELNAME_B=${MODEL%.json}
            FITFILE_B="$DIR/${MODELNAME_B}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME_B}.RData" -o $FITFILE_B;
            Rscript ../scripts/gather-results.R --fit $FITFILE_B --sim $DIR/sim.RData --json > ${FITFILE_B%RData}json;

            # BUT LETS SEE IF ANYTHING ELSE IS SIGNIFICANT
            MODELNAME_C="cpg-plus-1"
            FITFILE_C="$DIR/${MODELNAME_C}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL_C -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME_C}.RData" -o $FITFILE_C;
            Rscript ../scripts/gather-results.R --fit $FITFILE_C --sim $DIR/sim.RData --json > ${FITFILE_C%RData}json;

            # HOW BOUT WITH JUST ONE SPURIOUS ADDITION?
            MODELNAME_D="cpg-plus-2"
            FITFILE_D="$DIR/${MODELNAME_D}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL_D -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "genmatrices/genmatrix-${LONGWIN}-${MODELNAME_D}.RData" -o $FITFILE_D;
            Rscript ../scripts/gather-results.R --fit $FITFILE_D --sim $DIR/sim.RData --json > ${FITFILE_D%RData}json;

        done;
#     ) &
# done

../scripts/collect-params-results.R $DIR/*-fit*.json > cpg-plus-epsilon_results.tsv

wait;
