#!/bin/bash

set -eu

if [ $# -lt 4 ]
then
    echo "Usage: $0 (modlefile) (longwin) (shortwin) (leftwin)"
    exit 0
fi

MODELFILE=$1
LONGWIN=$2
SHORTWIN=$3
LEFTWIN=$4

MAXIT=500

DIRS="RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx"
TYPES="CTCF_binding_site  enhancer  open_chromatin_region  promoter_flanking_region"

# get the correct base frequencies in there:
MODELNAME=${MODELFILE%.json}

TMER="${LONGWIN}-${SHORTWIN}-${LEFTWIN}"

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

modelsetup $MODELFILE

collapse () {
    DIR=$1
    COLLAPSE=" library(contextual); f <- file.path(dir, \"9-5-2/total.counts.gz\"); g <- file.path(dirname(dirname(f)), \"${TMER}\", \"total.counts.gz\"); c <- read.counts(f, leftwin=1); d <- projectcounts(c, new.leftwin=${LEFTWIN}, new.shortwin=${SHORTWIN}, new.longwin=${LONGWIN}); dir.create(dirname(g), recursive=TRUE, showWarnings=FALSE); write.counts(d, file=g); "
    Rscript -e "dir=\"$DIR\"" -e "$COLLAPSE"
}

for d in $DIRS; do for t in $TYPES; do 
    if ! [ -e $d/$t/total.counts.gz ]
    then
        collapse $d/$t
    fi
done; done


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

# create genmatrix
GENMATFILE="genmatrices/${MODELNAME}-model-genmatrix-${LONGWIN}.RData"
if [ ! -e $GENMATFILE ]
then
    ../scripts/make-genmat.R -c models/$MODELFILE -w $LONGWIN -o $GENMATFILE
fi


# fit model
FITFILES=$(fitmodel $MODELFILE $TMER)

(for d in $DIRS; do for t in $TYPES; do if [ -e $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv ]; then
    echo $d $t; 
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    echo   "...."
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
fi; done; done)




