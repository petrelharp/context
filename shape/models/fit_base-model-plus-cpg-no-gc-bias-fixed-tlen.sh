#!/bin/bash

set -eu

DIRS="RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx"
TYPES="CTCF_binding_site  enhancer  open_chromatin_region  promoter_flanking_region"

## Base model:
# get the correct base frequencies in there:

modelsetup () {
    for d in $DIRS; do for t in $TYPES; do 
        CG=$(grep "$(basename $d)/$t" RegulatoryFeature-regions-from-axt/initfreqs.txt | cut -f 2 -d ' ')
        # JSON needs leading 0.'s ...
        A=0$(bc -l <<< "scale=4; (1-$CG)/2.0")
        C=0$(bc -l <<< "scale=4; $CG/2.0")
        G=0$(bc -l <<< "scale=4; $CG/2.0")
        T=0$(bc -l <<< "scale=4; (1-$CG)/2.0")
        cat models/$1 | sed -e 's_genmatrices_../../../genmatrices_' | sed -e "s/0.25, 0.25, 0.25, 0.25/$A, $C, $G, $T/" > $d/$t/$1; 
    done; done
}

MAXIT=500

# for fitting models, later
fitmodel () {
    MODELFILE=$1
    TMER=$2
    for d in $DIRS; do for t in $TYPES; do 
        ( FITFILE=$d/$t/$TMER/${MODELFILE%.json}-fit.RData
        if ! [ -e $FITFILE ]
        then
            echo "";
        fi
        ../scripts/fit-tree-model.R -i $d/$t/$TMER/total.counts.gz -c $d/$t/$MODELFILE -o $FITFILE -x $MAXIT
        ../scripts/gather-results.R --json -f $FITFILE > ${FITFILE%RData}json
        ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
        ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 
        ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 
        echo $FITFILE ) &
    done; done
    wait
}

## Worries that GC bias may introduce nonidentifiability
## ...  as above but no GC bias

# get the correct base frequencies in there:
MODELFILE="base-model-plus-cpg-no-gc-bias-fixed-tlen.json"
MODELNAME=${MODELFILE%.json}
modelsetup $MODELFILE

# create genmatrix
LONGWIN=4
GENMATFILE="genmatrices/${MODELNAME}-model-genmatrix-${LONGWIN}.RData"
if [ ! -e $GENMATFILE ]
then
    ../scripts/make-genmat.R -c models/$MODELFILE -w $LONGWIN -o $GENMATFILE
fi


# fit model
TMER="4-2-1"
FITFILES=$(fitmodel $MODELFILE $TMER)

(for d in $DIRS; do for t in $TYPES; do if [ -e $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv ]; then
    echo $d $t; 
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    echo   "...."
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
fi; done; done)



