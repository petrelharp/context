#!/bin/bash

set -eu

if [ $# -lt 4 ]
then
    echo "Run MCMC beginning from previous stopping point.\n Usage: $0 (modelfile) (longwin) (shortwin) (leftwin)"
    exit 0
fi

MODELFILE=$1
LONGWIN=$2
SHORTWIN=$3
LEFTWIN=$4
NBATCH=$5

DIRS="RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx"
TYPES="CTCF_binding_site  enhancer  open_chromatin_region  promoter_flanking_region"

# Find the previous config file
MODELNAME=$(basename ${MODELFILE%.json})

TMER="${LONGWIN}-${SHORTWIN}-${LEFTWIN}"

echo "Copying over models $MODELFILE to subdirectories..."
for d in $DIRS; do for t in $TYPES; do 
    echo "... $d/$t/$MODELFILE"
    CG=$(grep "$(basename $d)/$t" RegulatoryFeature-regions-from-axt/initfreqs.txt | cut -f 2 -d ' ')
    # JSON needs leading 0.'s ...
    A=0$(bc -l <<< "scale=4; (1-$CG)/2.0")
    C=0$(bc -l <<< "scale=4; $CG/2.0")
    G=0$(bc -l <<< "scale=4; $CG/2.0")
    T=0$(bc -l <<< "scale=4; (1-$CG)/2.0")
    cat $MODELFILE | sed -e 's_genmatrices_../../../genmatrices_' | sed -e "s/0.25, 0.25, 0.25, 0.25/$A, $C, $G, $T/" > $d/$t/${MODELNAME}.json; 
done; done


for d in $DIRS; do for t in $TYPES; do 

    FITFILE=$d/$t/$TMER/${MODELNAME}-fit.RData
    THIS_MODELFILE=$d/$t/${MODELNAME}.json
    if [ -e $FITFILE -a -e $THIS_MODELFILE ]
    then
        ( echo "MCMC on $FITFILE !";
          JOBID=$RANDOM
          OUTFILE=${FITFILE%.RData}-mcmc-${JOBID}.RData
          ../scripts/mcmc-tree-model.R -i $FITFILE -c $THIS_MODELFILE -b $NBATCH -j $JOBID -o $OUTFILE
          ../scripts/gather-results.R --json -f $OUTFILE > ${OUTFILE%RData}json
        echo $OUTFILE ) &
    else
        echo "Can't find files $FITFILE or $THIS_MODELFILE"
    fi
done; done
wait




