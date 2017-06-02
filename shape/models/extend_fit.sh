#!/bin/bash

set -eu

if [ $# -lt 4 ]
then
    echo "Run beginning from previous stopping point.\n Usage: $0 (modelfile) (longwin) (shortwin) (leftwin)"
    exit 0
fi

MODELFILE=$1
LONGWIN=$2
SHORTWIN=$3
LEFTWIN=$4

MAXIT=500

DIRS="RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx"
TYPES="CTCF_binding_site  enhancer  open_chromatin_region  promoter_flanking_region"

# Find the previous config file
MODELNAME=$(basename ${MODELFILE%.json})

TMER="${LONGWIN}-${SHORTWIN}-${LEFTWIN}"

for d in $DIRS; do for t in $TYPES; do 

    DATE=$(date +%d-%m-%Y-%H-%M-%S)
    FITFILE=$d/$t/$TMER/${MODELNAME}-fit.RData
    THIS_MODELFILE=$d/$t/${MODELNAME}.json
    TMER_MODELFILE=$d/$t/$TMER/${MODELNAME}.json
    echo "Updating $THIS_MODELFILE to $TMER_MODELFILE from $FITFILE"
    if [ -e $FITFILE -a -e $THIS_MODELFILE ]
    then
        ../scripts/update-config.R -m $FITFILE -c $THIS_MODELFILE -o $TMER_MODELFILE

        OLD_FITFILE=$d/$t/$TMER/${MODELNAME}-fit_${DATE}.RData
        echo "Moving $FITFILE to $OLD_FITFILE"
        cp -p $FITFILE $OLD_FITFILE

        ( echo "Fitting!";
          ../scripts/fit-tree-model.R -i $d/$t/$TMER/total.counts.gz -c $TMER_MODELFILE -o $FITFILE -x $MAXIT
          ../scripts/gather-results.R --json -f $FITFILE > ${FITFILE%RData}json
          ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
          ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 
          ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 
        echo $FITFILE ) &
    else
        echo "Can't find files $FITFILE or $THIS_MODELFILE"
    fi
done; done
wait

echo "The residuals:"
echo "--------------"
(for d in $DIRS; do for t in $TYPES; do if [ -e $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv ]; then
    echo $d $t; 
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    echo   "...."
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
fi; done; done)




