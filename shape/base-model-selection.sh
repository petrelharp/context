#!/bin/bash

DIRS="RegulatoryFeature-regions-from-axt/noOverlap-knownGeneTx RegulatoryFeature-regions-from-axt/overlap-knownGeneTx"
TYPES="CTCF_binding_site  enhancer  open_chromatin_region  promoter_flanking_region"

## Base model:
# get the correct base frequencies in there:
# NOTE: there is NO difference in base frequencies between the two taxa

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

modelsetup base-model.json

# for fitting models, later
MAXIT=500

fitmodel () {
    MODELFILE=$1
    TMER=$2
    for d in $DIRS; do for t in $TYPES; do 
        FITFILE=$d/$t/$TMER/${MODELFILE%.json}-fit.RData
        if ! [ -e $FITFILE ]
        then
            ../scripts/fit-tree-model.R -i $d/$t/$TMER/total.counts.gz -c $d/$t/$MODELFILE -o $FITFILE -x $MAXIT
        fi
        ../scripts/gather-results.R --json -f $FITFILE > ${FITFILE%RData}json
        ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l1.tsv -w 2 -s 1 -l 1 --pretty 
        ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.2.1.l0.tsv -w 2 -s 1 -l 0 --pretty 
        ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.3.1.l1.tsv -w 3 -s 1 -l 1 --pretty 
        echo $FITFILE
    done; done
}

getresids () {
    LONGWIN=$1
    SHORTWIN=$2
    LEFTWIN=$3
    for d in $DIRS; do for t in $TYPES; do 
        for FITFILE in $d/$t/*/*-fit.RData;
        do
            OUTFILE=${FITFILE%.RData}-resids.${LONGWIN}.${SHORTWIN}.l${LEFTWIN}.tsv
            if ! [ -e $OUTFILE ]
            then
                if [ $# -gt 3 ]
                then
                    # longer Tmers
                    COUNTTMER=$4
                    COUNTFILE=$d/$t/${COUNTTMER}/total.counts.gz
                    FITFILENAME=$(basename $FITFILE)
                    MODELFILE=$(dirname $(dirname $FITFILE))/${FITFILENAME%-fit.RData}.json
                    COUNTARG="--countfile $COUNTFILE --config $MODELFILE"
                else
                    COUNTARG=""
                fi
                ../scripts/compute-resids.R -i $FITFILE -o ${FITFILE%.RData}-resids.${LONGWIN}.${SHORTWIN}.l${LEFTWIN}.tsv -w ${LONGWIN} -s ${SHORTWIN} -l ${LEFTWIN} --pretty $COUNTARG
                echo $OUTFILE
            fi
        done
    done; done
}

# collapse down to 6-3-2
collapse () {
    DIR=$1
    # note leftwin must match input counts here
    COLLAPSE=" library(contextual); f <- file.path(dir, \"9-5-2/total.counts.gz\"); g <- file.path(dirname(dirname(f)), \"6-3-2\", \"total.counts.gz\"); c <- read.counts(f, leftwin=2); d <- projectcounts(c, new.leftwin=2, new.shortwin=3, new.longwin=6); dir.create(dirname(g), recursive=TRUE, showWarnings=FALSE); write.counts(d, file=g); "
    Rscript -e "dir=\"$DIR\"" -e "$COLLAPSE"
}

for d in $DIRS; do for t in $TYPES; do 
    collapse $d/$t
done; done

# create genmatrix
LONGWIN=5
for MODELFILE in models/*json
do
    MODELNAME=$(basename ${MODELFILE%.json})
    if [ ! -e $GENMATFILE ]
    then
        ../scripts/make-genmat.R -c $MODELFILE -w $LONGWIN -o genmatrices/
    fi
done


# fit model
TMER="6-3-2"
for d in $DIRS; do for t in $TYPES; do 
    ../scripts/fit-tree-model.R -i $d/$t/$TMER/total.counts.gz -c $d/$t/base-model.json -o $d/$t/$TMER/base-model-fit.RData 
done; done

######## SHORTER

# collapse down to 4-2-1
collapsemore () {
    DIR=$1
    COLLAPSE=" library(contextual); f <- file.path(dir, \"7-5-1/total.counts.gz\"); g <- file.path(dirname(dirname(f)), \"4-2-1\", \"total.counts.gz\"); c <- read.counts(f, leftwin=1); d <- projectcounts(c, new.leftwin=1, new.shortwin=2, new.longwin=4); dir.create(dirname(g), recursive=TRUE, showWarnings=FALSE); write.counts(d, file=g); "
    Rscript -e "dir=\"$DIR\"" -e "$COLLAPSE"
}

for d in $DIRS; do for t in $TYPES; do 
    collapsemore $d/$t
done; done

# create genmatrix
LONGWIN=4
../scripts/make-genmat.R -c base-model.json -w $LONGWIN -o genmatrices/base-model-genmatrix-${LONGWIN}.RData

# fit model

TMER="4-2-1"
fitmodel base-model.json $TMER

# clearly, need to add CpG
(for d in $DIRS; do for t in $TYPES; do 
    echo $d $t; 
    head -n 5 $d/$t/$TMER/base-model-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    tail -n 5 $d/$t/$TMER/base-model-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    echo   "...."
    head -n 5 $d/$t/$TMER/base-model-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
    tail -n 5 $d/$t/$TMER/base-model-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
done; done)

# ... and that's the main thing
(for d in $DIRS; do for t in $TYPES; do 
    echo $d $t; 
    tail -n 5 $d/$t/$TMER/base-model-fit-resids.3.1.l1.tsv | awk '{print $2"\t."$3".\t"$7}';
done; done)


######################################
## Revised model with CG->*
## including making base model strand-symmetric and setting rates to near inferred (on first pass)

# get the correct base frequencies in there:
MODELFILE="base-model-plus-cpg.json"
MODELNAME=${MODELFILE%.json}
modelsetup $MODELFILE

# create genmatrix
LONGWIN=4
../scripts/make-genmat.R -c models/$MODELFILE -w $LONGWIN -o genmatrices/${MODELNAME}-model-genmatrix-${LONGWIN}.RData

# fit model
TMER="4-2-1"
fitmodel $MODELFILE $TMER

(for d in $DIRS; do for t in $TYPES; do if [ -e $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv ]; then
    echo $d $t; 
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l0.tsv | awk '{print $2"\t."$3"\t"$7}';
    echo   "...."
    head -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
    tail -n 5 $d/$t/$TMER/${MODELNAME}-fit-resids.2.1.l1.tsv | awk '{print $2"\t"$3".\t"$7}';
fi; done; done)


######################################
## Worries that GC bias may introduce nonidentifiability
## ...  as above but no GC bias

# get the correct base frequencies in there:
MODELFILE="base-model-plus-cpg-no-gc-bias.json"
MODELNAME=${MODELFILE%.json}
modelsetup $MODELFILE

# create genmatrix
LONGWIN=4
../scripts/make-genmat.R -c models/$MODELFILE -w $LONGWIN -o genmatrices/${MODELNAME}-model-genmatrix-${LONGWIN}.RData

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

(for d in $DIRS; do for t in $TYPES; do if [ -e $d/$t/$TMER/${MODELNAME}-fit-resids.3.1.l1.tsv ]; then
    echo $d $t; 
    head -n 10 $d/$t/$TMER/${MODELNAME}-fit-resids.3.1.l1.tsv | awk '{print $2"\t."$3".\t"$7}';
    tail -n 10 $d/$t/$TMER/${MODELNAME}-fit-resids.3.1.l1.tsv | awk '{print $2"\t."$3".\t"$7}';
fi; done; done)




##########################
# Turns out that things run off into a situation of having the 'derived' branch
# much longer than the 'reference' branch, with mutation rates going down to
# compensate.  Strangely, this has much better likelihood. Introduced the
# "-fixed-tlens' models to account for this, especially as in cases where this
# didn't happen, branches were almost identical.

## GC-bias is not substantially bettering the likelihood.


## Moved stuff over to fit_model.sh:

../scripts/cluster/submit-8-jobs.sh models/fit_model.sh cpg-plus-epsilon-separate-branches.json 6 3 2
../scripts/cluster/submit-8-jobs.sh models/fit_model.sh cpg-plus-epsilon-separate-branches.json 7 5 1
../scripts/cluster/submit-8-jobs.sh models/fit_model.sh shape-model-CpG-only-MGW.json 7 5 1
../scripts/cluster/submit-8-jobs.sh models/fit_model.sh shape-model-CpG-only-random-values.json 7 5 1

../scripts/cluster/submit-8-jobs.sh models/fit_model.sh models/base-model-plus-cpg-and-CPD.json 4 2 1
../scripts/cluster/submit-8-jobs.sh models/fit_model.sh models/base-model-plus-cpg-and-CPD.json 6 3 2
../scripts/cluster/submit-8-jobs.sh models/fit_model.sh models/base-model-plus-cpg-and-CPD.json 7 3 2


#################################
### collect results

# look at all the results
../scripts/collect-params-results.R RegulatoryFeature-regions-from-axt/**/*fit.json > model-selection-results.tsv

paste <(cat model-selection-results.tsv | cut -f 1 | sed -e 's_.*/__') <( cat model-selection-results.tsv | cut -f 4,21,22-25) | column -t

# look at all the bottom/top 1-1-0 resids:
(for x in RegulatoryFeature-regions-from-axt/**/*resids.1.1.l0.tsv; do y=$(basename $x | sed -e 's/-resids.*//'); z=$(head -n 2 $x | tail -n 1); echo $y $z; done) | column -t
(for x in RegulatoryFeature-regions-from-axt/**/*resids.1.1.l0.tsv; do y=$(basename $x | sed -e 's/-resids.*//'); z=$(tail -n 1 $x); echo $y $z; done) | column -t

# look at all the top 4-2-1 resids:
(for x in RegulatoryFeature-regions-from-axt/**/*resids.4.2.l1.tsv; do y=$(basename $x | sed -e 's/-resids.*//'); z=$(head -n 2 $x | tail -n 3); echo $y $z; done) | column -t
(for x in RegulatoryFeature-regions-from-axt/**/*resids.4.2.l1.tsv; do y=$(basename $x | sed -e 's/-resids.*//'); z=$(tail -n 3 $x); echo $y $z; done) | column -t


showresids () {
    LONGWIN=$1; SHORTWIN=$2; LEFTWIN=$3; NPATS=${4:-1}
    # look at all the top 4-2-1 resids:
    (for x in RegulatoryFeature-regions-from-axt/**/*resids.${LONGWIN}.${SHORTWIN}.l${LEFTWIN}.tsv; do y=$(basename $x | sed -e 's/-resids.*//'); z=$(head -n $((NPATS+1)) $x | tail -n +2); echo $y $z; done) | column -t
    (for x in RegulatoryFeature-regions-from-axt/**/*resids.${LONGWIN}.${SHORTWIN}.l${LEFTWIN}.tsv; do y=$(basename $x | sed -e 's/-resids.*//'); z=$(tail -n $NPATS $x); echo $y $z; done) | column -t
}
