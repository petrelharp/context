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

# where to put stuff
BASEDIR=cpg-plus-epsilon
mkdir -p $BASEDIR

TRUE_MODEL="$BASEDIR/cpg-plus-epsilon.json"
# single-base-only model initialized for mean sequence divergence
BASEMODEL="$BASEDIR/base-model.json"
# single-base-only with CpG, initialized to the previous model
MODEL_A="$BASEDIR/cpg-only.json"
# as above but with Harris's mutations added, initialized to previous model
# note this is structurally the same as TRUE_MODEL
MODEL_B="$BASEDIR/cpg-and-harris.json"

## --------

MODEL2="$BASEDIR/base-plus-2mers.json"
MODEL3="$BASEDIR/base-plus-3mers.json"

MODEL_C="cpg-plus-epsilon/cpg-plus-1.json"
MODEL_D="cpg-plus-epsilon/cpg-plus-2.json"

SUBDIR="simseqs/cpg-plus"
TLEN=1.0

# precompute generator matrices:
GMDIR=$BASEDIR/genmatrices
mkdir -p $GMDIR
for THISMODEL in $TRUE_MODEL $BASEMODEL $MODEL_A $MODEL_B
do
    MODELNAME=$(basename ${THISMODEL%.json})
    for LONGWIN in 4 5
    do
        GENMAT="$GMDIR/genmatrix-${LONGWIN}-${MODELNAME}.RData"
        if [[ ! -e $GENMAT ]]
        then
            Rscript ../scripts/make-genmat.R -c $THISMODEL -w ${LONGWIN} -o $GENMAT
        fi
    done
done

SEQLEN=1000000


# for N in $(seq $NRUNS)
# do
#     (
#        DIR=$BASEDIR/$(printf "%05g" $RANDOM);

fitmodel ()
{
    THISMODEL=$1
    MODELNAME=$(basename ${THISMODEL%.json})
    FITFILE="$DIR/${MODELNAME}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
    LOGFILE=${FITFILE%RData}log
    Rscript ../scripts/fit-model.R -c $THISMODEL -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "$GMDIR/genmatrix-${LONGWIN}-${MODELNAME}.RData" -o $FITFILE &>$LOGFILE
    Rscript ../scripts/gather-results.R --fit $FITFILE --sim $DIR/sim.RData --json > ${FITFILE%RData}json 2>$LOGFILE
    echo $FITFILE
}

DIR=$BASEDIR/32140

        echo "Simulation $N, in $DIR";
        mkdir -p "$DIR";
        # simulate up some sequence for testing;
        if [ ! -f $DIR/sim.RData ]
        then
            Rscript ../scripts/sim-seq.R -c $TRUE_MODEL -t $TLEN -s $SEQLEN -d "$DIR" -o "sim.RData";
        fi

        # for LONGWIN in 4 5;
        # do
        LONGWIN=5
            # this gets 1-2-1 and 1-3-1
            SHORTWIN=$(( (LONGWIN+1)/2 ));
            LEFTWIN=$(( (LONGWIN-SHORTWIN)/2 ));

            # count the Tmers;
            Rscript ../scripts/count-seq.R -i $DIR/sim.RData -w $LONGWIN -s $SHORTWIN -l $LEFTWIN;

            # fit the basic model;
            BASIC_FIT=$(fitmodel $BASEMODEL)

            # compute residuals
            ../scripts/compute-resids.R -i $BASIC_FIT -o ${BASIC_FIT%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
            ../scripts/compute-resids.R -i $BASIC_FIT -o ${BASIC_FIT%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 

            # add in CpG - initial parameters copied from previous fit
            FIT_A=$(fitmodel $MODEL_A)

            # compute residuals
            ../scripts/compute-resids.R -i $FIT_A -o ${FIT_A%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
            ../scripts/compute-resids.R -i $FIT_A -o ${FIT_A%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 
            ../scripts/compute-resids.R -i $FIT_A -o ${FIT_A%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 

            ## -------------

            # add in Harris - initial parameters copied from previous fit
            FIT_B=$(fitmodel $MODEL_B)

            # compute residuals
            ../scripts/compute-resids.R -i $FIT_B -o ${FIT_B%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
            ../scripts/compute-resids.R -i $FIT_B -o ${FIT_B%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 
            ../scripts/compute-resids.R -i $FIT_B -o ${FIT_B%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 


            # now with all 2mers
            MODELNAME2="${MODEL2%.json}"
            FITFILE2="$DIR/2mer-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL2 -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "$GMDIR/genmatrix-${LONGWIN}-${MODELNAME2}.RData" -o $FITFILE2;
            Rscript ../scripts/gather-results.R --fit $FITFILE2 --sim $DIR/sim.RData --json > ${FITFILE2%RData}json;

            # CLEAR FROM RESIDUALS: now with CpG
            MODELNAME_A=${MODEL_A%.json}
            FITFILE_A="$DIR/${MODELNAME_A}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL_A -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "$GMDIR/genmatrix-${LONGWIN}-${MODELNAME_A}.RData" -o $FITFILE_A;
            Rscript ../scripts/gather-results.R --fit $FITFILE_A --sim $DIR/sim.RData --json > ${FITFILE_A%RData}json;
            ../scripts/compute-resids.R -i $FITFILE_A -o ${FITFILE_A%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 

            # AGAIN CLEAR FROM RESIDUALS
            MODELNAME_B=${MODEL%.json}
            FITFILE_B="$DIR/${MODELNAME_B}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "$GMDIR/genmatrix-${LONGWIN}-${MODELNAME_B}.RData" -o $FITFILE_B;
            Rscript ../scripts/gather-results.R --fit $FITFILE_B --sim $DIR/sim.RData --json > ${FITFILE_B%RData}json;
            ../scripts/compute-resids.R -i $FITFILE_B -o ${FITFILE_B%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 

            # BUT LETS SEE IF ANYTHING ELSE IS SIGNIFICANT
            MODELNAME_C="cpg-plus-1"
            FITFILE_C="$DIR/${MODELNAME_C}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL_C -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "$GMDIR/genmatrix-${LONGWIN}-${MODELNAME_C}.RData" -o $FITFILE_C;
            Rscript ../scripts/gather-results.R --fit $FITFILE_C --sim $DIR/sim.RData --json > ${FITFILE_C%RData}json;

            # HOW BOUT WITH JUST ONE SPURIOUS ADDITION?
            MODELNAME_D="cpg-plus-2"
            FITFILE_D="$DIR/${MODELNAME_D}-fit-${LONGWIN}-${SHORTWIN}-l${LEFTWIN}.RData";
            Rscript ../scripts/fit-model.R -c $MODEL_D -t $TLEN -i "$DIR/sim-${LONGWIN}-root-${SHORTWIN}-tip-l${LEFTWIN}-shift0.counts" -m "$GMDIR/genmatrix-${LONGWIN}-${MODELNAME_D}.RData" -o $FITFILE_D;
            Rscript ../scripts/gather-results.R --fit $FITFILE_D --sim $DIR/sim.RData --json > ${FITFILE_D%RData}json;

        # done;
#     ) &
# done

../scripts/collect-params-results.R $DIR/*-fit*.json > cpg-plus-epsilon_results.tsv

wait;
