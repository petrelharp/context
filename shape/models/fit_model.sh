#!/bin/bash

set -eu

if [ $# -lt 4 ]
then
    echo "Usage: $0 (modelfile) (longwin) (shortwin) (leftwin)"
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
MODELNAME=$(basename ${MODELFILE%.json})

TMER="${LONGWIN}-${SHORTWIN}-${LEFTWIN}"

# get the correct base frequencies in there:

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

collapse () {
    DIR=$1
    # note that leftwin must match input counts here
    COLLAPSE="library(contextual); f <- file.path(dir, \"9-5-2/total.counts.gz\"); g <- file.path(dirname(dirname(f)), \"${TMER}\", \"total.counts.gz\"); c <- read.counts(f, leftwin=2); d <- projectcounts(c, new.leftwin=${LEFTWIN}, new.shortwin=${SHORTWIN}, new.longwin=${LONGWIN}); dir.create(dirname(g), recursive=TRUE, showWarnings=FALSE); write.counts(d, file=g); "
    Rscript -e "dir=\"$DIR\"" -e "$COLLAPSE"
}

echo "Building count files."
for d in $DIRS; do for t in $TYPES; do 
    if ! [ -e $d/$t/total.counts.gz ]
    then
        echo " ... $d/$t/total.counts.gz"
        collapse $d/$t &  # in parallel
    else
        echo " ... already exists: $d/$t/total.counts.gz"
    fi
done; done

wait

# create genmatrix
GENMATFILE="genmatrices/${MODELNAME}-model-genmatrix-${LONGWIN}.RData"
if [ ! -e $GENMATFILE ]
then
    echo "Making generator matrix $GENMATFILE ..."
    ../scripts/make-genmat.R -c $MODELFILE -w $LONGWIN -o $GENMATFILE
fi

# for fitting models, later
fitmodel () {
    MODELNAME=$1
    TMER=$2
    for d in $DIRS; do for t in $TYPES; do 
        ( FITFILE=$d/$t/$TMER/${MODELNAME}-fit.RData
        echo "... $FITFILE"
        if ! [ -e $FITFILE ]
        then
            echo "Fitting!";
            ../scripts/fit-tree-model.R -i $d/$t/$TMER/total.counts.gz -c $d/$t/${MODELNAME}.json -o $FITFILE -x $MAXIT
            ../scripts/gather-results.R --json -f $FITFILE > ${FITFILE%RData}json
            ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
            ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 
            ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 
        else
            echo "Fit already exists!";
        fi
        echo $FITFILE ) &
    done; done
    wait
}

# fit model
echo "Fitting models to $TMER ..."
FITFILES=$(fitmodel $MODELNAME $TMER)
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




